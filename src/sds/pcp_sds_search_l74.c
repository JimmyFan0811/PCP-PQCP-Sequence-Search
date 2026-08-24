#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define L 74

// --------- 24 個 orbit（每個大小 3），不含 0 與 37 ---------
#define NUM_ORBITS 24
#define ORBIT_SIZE 3

// 每個 orbit 的三個 index
int orbits[NUM_ORBITS][ORBIT_SIZE] = {
    {1, 47, 63},
    {3, 41, 67},
    {2, 20, 52},
    {4, 30, 40},
    {5, 13, 19},
    {6,  8, 60},
    {7, 33, 71},
    {9, 49, 53},
    {10, 26, 38},
    {11, 27, 73},
    {12, 16, 46},
    {14, 66, 68},
    {15, 39, 57},
    {17, 35, 59},
    {18, 24, 32},
    {21, 25, 65},
    {22, 54, 72},
    {23, 43, 45},
    {28, 58, 62},
    {29, 31, 51},
    {34, 44, 70},
    {36, 48, 64},
    {42, 50, 56},
    {55, 61, 69}
};

// --------- 參數 ---------
#define BATCH_A 5000    // 一批最多 1000 個 A
const int KA = 12;           // A: 選 12 個 orbit 為 -1
const int KB = 10;           // B: 選 10 個 orbit 為 -1
const unsigned long long TOTAL_A_COMBOS = 2704156ULL; // C(24,12)

// --------- 輔助函式：PACF / 印出 / 組合 ---------

int periodic_autocorrelation(int *seq, int shift) {
    int sum = 0;
    for (int i = 0; i < L; i++) {
        int j = i + shift;
        if (j >= L) j -= L;
        sum += (seq[i] == seq[j]) ? 1 : -1;
    }
    return sum;
}

void print_sequence(int *seq, const char *name) {
    printf("%s = { ", name);
    for (int i = 0; i < L; i++) {
        putchar(seq[i] == 1 ? '+' : '-');
    }
    printf(" };\n");
}

void write_seq_to_file(FILE *fp, int *seq, const char *name) {
    fprintf(fp, "%s = { ", name);
    for (int i = 0; i < L; i++) {
        fputc(seq[i] == 1 ? '+' : '-', fp);
    }
    fprintf(fp, " };\n");
}

// n 個物件中選 k 個的組合：idx[0..k-1] 維持遞增
// 初始化後，呼叫 next_combination 取得下一個組合
int next_combination(int *idx, int n, int k) {
    int i = k - 1;
    while (i >= 0 && idx[i] == n - k + i) {
        i--;
    }
    if (i < 0) return 0; // 沒有下一個組合

    idx[i]++;
    for (int j = i + 1; j < k; j++) {
        idx[j] = idx[j - 1] + 1;
    }
    return 1;
}

// 用 orbit 組合 combA 產生一個 A 序列
void build_A_from_comb(int *combA, int *seq_a) {
    // 全部先設成 +1
    for (int i = 0; i < L; i++) seq_a[i] = 1;
    // 0, 37 固定 +1
    seq_a[0]  = 1;
    seq_a[37] = 1;

    // 將 12 個 orbit 的位置設成 -1
    for (int t = 0; t < KA; t++) {
        int o = combA[t];
        for (int j = 0; j < ORBIT_SIZE; j++) {
            int pos = orbits[o][j];
            seq_a[pos] = -1;
        }
    }
}

// 檢查指定的 A (pacf_a, max_abs_a) 與 B (pacf_b) 是否為 PCP
int verify_pcp_with_precomputed(const int *pacf_a, int max_abs_a, const int *pacf_b) {
    // 先用 |PACF_B| <= max_abs_a 做 pruning
    for (int u = 1; u <= L / 2; u++) {
        if (abs(pacf_b[u]) > max_abs_a) {
            return 0;
        }
    }
    // 完整 PCP 條件：PACF_A(u) + PACF_B(u) = 0
    for (int u = 1; u <= L / 2; u++) {
        if (pacf_a[u] + pacf_b[u] != 0) {
            return 0;
        }
    }
    return 1;
}

int main(void) {
    // 一批 A 的資料
    static int seqA_batch[BATCH_A][L];
    static int pacfA_batch[BATCH_A][L];
    static int maxAbsA_batch[BATCH_A];

    // B 序列與 PACF(B)
    int seq_b[L];
    int pacf_b[L];

    unsigned long long checked_pairs = 0ULL;  // A-B 配對數量
    unsigned long long found_pairs   = 0ULL;  // 找到幾組 PCP
    unsigned long long A_done        = 0ULL;  // 已處理的 A 組數

    FILE *fp = fopen("pcp_pairs_len74_speed1.txt", "w");
    if (!fp) {
        perror("cannot open pcp_pairs_len74.txt");
        return 1;
    }

    printf("Start searching PCPs for binary length L = %d\n", L);
    printf("Total A combinations = %llu (C(24,12))\n", TOTAL_A_COMBOS);
    printf("Batch size for A = %d\n", BATCH_A);
    printf("A: choose %d orbits (out of %d) as -1, plus {0,37}=+1.\n", KA, NUM_ORBITS);
    printf("B: choose %d orbits as -1, plus {0} or {37} as extra -1.\n", KB);

    // 初始化第一個 A 組合 0..11
    int combA[KA];
    for (int i = 0; i < KA; i++) combA[i] = i;

    // ======= 外圈：不斷產生新的 A 批次 =======
    while (A_done < TOTAL_A_COMBOS) {
        int batch_count = 0;  // 這一批實際有幾個 A

        // ---- STEP 1: 產生這一批的 A、PACF(A) ----
        while (batch_count < BATCH_A && A_done < TOTAL_A_COMBOS) {
            int *seq_a   = seqA_batch[batch_count];
            int *pacf_a  = pacfA_batch[batch_count];

            build_A_from_comb(combA, seq_a);

            // sanity check: A 的和應該是 2 (38 個 +1，36 個 -1)
            int sumA = 0;
            for (int i = 0; i < L; i++) sumA += seq_a[i];
            if (sumA != 2) {
                printf("Warning: A sum != 2 (sumA=%d)\n", sumA);
            }

            // 計算 A 的 PACF 與 max |PACF|
            int max_abs = 0;
            for (int u = 0; u < L; u++) {
                pacf_a[u] = periodic_autocorrelation(seq_a, u);
                if (u > 0 && abs(pacf_a[u]) > max_abs) {
                    max_abs = abs(pacf_a[u]);
                }
            }
            maxAbsA_batch[batch_count] = max_abs;

            batch_count++;
            A_done++;

            // 準備下一個 A 的組合
            if (!next_combination(combA, NUM_ORBITS, KA)) {
                // 沒有更多 A 了
                break;
            }
        }

        if (batch_count == 0) break;

        // ---- STEP 2: 針對這一批 A，枚舉所有 B ----
        int combB[KB];
        for (int i = 0; i < KB; i++) combB[i] = i;
        int has_more_B = 1;

        while (has_more_B) {
            // extra_singleton = 0 → index 0 為 -1
            // extra_singleton = 1 → index 37 為 -1
            for (int extra_singleton = 0; extra_singleton < 2; extra_singleton++) {

                // 建構 B 序列：全部先 +1
                for (int i = 0; i < L; i++) seq_b[i] = 1;

                // 10 個 orbit 設成 -1
                for (int t = 0; t < KB; t++) {
                    int o = combB[t];
                    for (int j = 0; j < ORBIT_SIZE; j++) {
                        int pos = orbits[o][j];
                        seq_b[pos] = -1;
                    }
                }

                // 單點 {0} 或 {37} 設成 -1
                if (extra_singleton == 0) {
                    seq_b[0] = -1;
                } else {
                    seq_b[37] = -1;
                }

                // sanity check: B 的和應該是 12 (43 個 +1，31 個 -1)
                int sumB = 0;
                for (int i = 0; i < L; i++) sumB += seq_b[i];
                if (sumB != 12) {
                    printf("Warning: B sum != 12 (sumB=%d)\n", sumB);
                }

                // 計算這個 B 的 PACF，一次即可給整批 A 用
                for (int u = 0; u < L; u++) {
                    pacf_b[u] = periodic_autocorrelation(seq_b, u);
                }

                // ---- STEP 3：對於這一批的所有 A 做檢查（OpenMP 平行化）----
                #pragma omp parallel for default(shared)
                for (int a_idx = 0; a_idx < batch_count; a_idx++) {
                    int is_pcp = verify_pcp_with_precomputed(
                        pacfA_batch[a_idx],
                        maxAbsA_batch[a_idx],
                        pacf_b
                    );

                    unsigned long long local_checked;
                    // checked_pairs++，並擷取目前的值
                    #pragma omp atomic capture
                    {
                        checked_pairs++;
                        local_checked = checked_pairs;
                    }

                    // 依舊每 1000000 次印一次進度（多執行緒下略有非同步，但整體 OK）
                    if (local_checked % 1000000ULL == 0ULL) {
                        #pragma omp critical
                        {
                            double percent = (double)A_done * 100.0 /
                                             (double)TOTAL_A_COMBOS;
                            printf("Checked %llu pairs so far (A combos: %llu / %llu, %.4f%%).\n",
                                   local_checked, A_done, TOTAL_A_COMBOS, percent);
                        }
                    }

                    if (is_pcp) {
                        unsigned long long this_pair_idx;
                        #pragma omp critical
                        {
                            found_pairs++;
                            this_pair_idx = found_pairs;

                            printf("\n=== Found PCP pair #%llu ===\n", this_pair_idx);
                            print_sequence(seqA_batch[a_idx], "A");
                            print_sequence(seq_b, "B");

                            fprintf(fp, "=== Found PCP pair #%llu ===\n", this_pair_idx);
                            write_seq_to_file(fp, seqA_batch[a_idx], "A");
                            write_seq_to_file(fp, seq_b, "B");
                            fprintf(fp, "\n");
                            fflush(fp);
                        }
                    }
                } // end parallel for over A
            }

            has_more_B = next_combination(combB, NUM_ORBITS, KB);
        }

        // 這一批 A 完成
        double percent_done = (double)A_done * 100.0 / (double)TOTAL_A_COMBOS;
        printf("\n--- Finished an A batch (current A combos: %llu / %llu, %.4f%%) ---\n",
               A_done, TOTAL_A_COMBOS, percent_done);
    }

    printf("\nSearch finished.\n");
    printf("Total A combinations processed = %llu (out of %llu)\n",
           A_done, TOTAL_A_COMBOS);
    printf("Total pairs checked   = %llu\n", checked_pairs);
    printf("Total PCP pairs found = %llu\n", found_pairs);

    fclose(fp);
    return 0;
}
