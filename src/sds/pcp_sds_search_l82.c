#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 82

// ============================
// Orbit definitions (order = 6)
// ============================

// 16 orbits of size 5
#define NUM_ORBITS5 16
int orbits5[NUM_ORBITS5][5] = {
    {1, 37, 57, 59, 51},
    {2, 74, 32, 36, 20},
    {3, 29, 7,  13, 71},
    {4, 66, 64, 72, 40},
    {5, 21, 39, 9, 49},
    {6, 58, 14, 60, 26},
    {8, 50, 46, 62, 80},
    {10, 42, 78, 16, 18},
    {11, 79, 53, 75, 69},
    {12, 28, 34, 52, 38},
    {15, 63, 35, 27, 65},
    {17, 55, 67, 19, 47},
    {22, 24, 76, 68, 56},
    {23, 31, 81, 45, 25},
    {30, 70, 48, 44, 54},
    {33, 73, 77, 61, 43}
};

// 2 trivial orbits of size 1
#define NUM_ORBITS1 2
int orbits1[NUM_ORBITS1] = {
    0, 41
};

// ============================
// Search parameters
// ============================

// A : 45 negatives  => 5*a5 + a1 = 45 => a5 = 9, a1 = 0
#define A_K5 9
#define A_K1 0

// B : 36 negatives  => 5*b5 + b1 = 36 => b5 = 7, b1 = 1
#define B_K5 7
#define B_K1 1

// Total A combos = C(16,9) = 11440
#define TOTAL_A_COMBOS 11440ULL

// Batch size for A
#define BATCH_A 1000ULL

// ============================
// Combination iterator
// ============================

int next_combination(int *idx, int n, int k) {
    int i = k - 1;
    while (i >= 0 && idx[i] == n - k + i) i--;
    if (i < 0) return 0;
    idx[i]++;
    for (int j = i + 1; j < k; j++)
        idx[j] = idx[j - 1] + 1;
    return 1;
}

// ============================
// PACF
// ============================

int periodic_autocorrelation(int *seq, int shift) {
    int sum = 0;
    for (int i = 0; i < L; i++) {
        int j = i + shift;
        if (j >= L) j -= L;
        sum += (seq[i] == seq[j]) ? 1 : -1;
    }
    return sum;
}

// ============================
// Print helpers
// ============================

void print_sequence(int *s, const char *name) {
    printf("%s = { ", name);
    for (int i = 0; i < L; i++)
        putchar(s[i] == 1 ? '+' : '-');
    printf(" };\n");
}

void write_sequence(FILE *fp, int *s, const char *name) {
    fprintf(fp, "%s = { ", name);
    for (int i = 0; i < L; i++)
        fputc(s[i] == 1 ? '+' : '-', fp);
    fprintf(fp, " };\n");
}

// ============================
// PCP check
// ============================

int verify_pcp_with_precomputed(const int *pacf_a, int max_abs_a,
                                const int *pacf_b) {
    for (int u = 1; u <= L / 2; u++) {
        if (abs(pacf_b[u]) > max_abs_a) return 0;
    }
    for (int u = 1; u <= L / 2; u++) {
        if (pacf_a[u] + pacf_b[u] != 0) return 0;
    }
    return 1;
}

// ============================
// main
// ============================

int main(void) {
    static int seqA_batch[BATCH_A][L];
    static int pacfA_batch[BATCH_A][L];
    static int maxAbsA_batch[BATCH_A];

    int seq_b[L];
    int pacf_b[L];

    unsigned long long checked_pairs = 0ULL;
    unsigned long long found_pairs   = 0ULL;
    unsigned long long A_done        = 0ULL;

    FILE *fp = fopen("pcp_pairs_len82_order6.txt", "w");
    if (!fp) {
        perror("cannot open output file");
        return 1;
    }

    printf("Start searching PCPs for L = %d (order = 6)\n", L);
    printf("Total A combinations = %llu\n", TOTAL_A_COMBOS);
    printf("Batch size = %llu\n", (unsigned long long)BATCH_A);

    // -----------------------------
    // A 的 combination 初始化
    // -----------------------------
    int A_comb5[A_K5];
    int A_has_next = 1;
    // 第一組: {0,1,...,A_K5-1}
    for (int i = 0; i < A_K5; i++) A_comb5[i] = i;

    while (A_done < TOTAL_A_COMBOS && A_has_next) {
        int batch_count = 0;

        // ========= 產生一批 A =========
        while (batch_count < (int)BATCH_A &&
               A_done < TOTAL_A_COMBOS &&
               A_has_next) {

            int *seq_a  = seqA_batch[batch_count];
            int *pacf_a = pacfA_batch[batch_count];

            // 先全部設為 +1
            for (int i = 0; i < L; i++) seq_a[i] = 1;

            // 將選到的 5-element orbits 設為 -1
            for (int t = 0; t < A_K5; t++) {
                int o = A_comb5[t];
                for (int k = 0; k < 5; k++)
                    seq_a[orbits5[o][k]] = -1;
            }
            // A_K1 = 0，這裡不用處理 trivial

            // 檢查 sum(A) 是否正確：82 - 2*45 = -8
            int sumA = 0;
            for (int i = 0; i < L; i++) sumA += seq_a[i];
            if (sumA != -8) {
                printf("Warning: A sum != -8 (sumA=%d)\n", sumA);
            }

            // 預算 PACF(A) 以及 max |PACF|
            int max_abs = 0;
            for (int u = 0; u < L; u++) {
                pacf_a[u] = periodic_autocorrelation(seq_a, u);
                if (u > 0 && abs(pacf_a[u]) > max_abs)
                    max_abs = abs(pacf_a[u]);
            }
            maxAbsA_batch[batch_count] = max_abs;

            batch_count++;
            A_done++;

            // 下一組 A 的 orbit 組合
            if (!next_combination(A_comb5, NUM_ORBITS5, A_K5)) {
                A_has_next = 0;
            }
        }

        if (batch_count == 0) break;

        // ========= 這一批 A，掃過所有 B =========
        // B: b5 = 7, b1 = 1
        int B_comb5[B_K5];

        for (int b1_idx = 0; b1_idx < NUM_ORBITS1; b1_idx++) {

            // 初始化 7-of-16 組合
            for (int i = 0; i < B_K5; i++) B_comb5[i] = i;
            int B_has_next = 1;

            while (B_has_next) {
                // 建立 B 序列
                for (int i = 0; i < L; i++) seq_b[i] = 1;

                // 非平凡軌道
                for (int t = 0; t < B_K5; t++) {
                    int o = B_comb5[t];
                    for (int k = 0; k < 5; k++)
                        seq_b[orbits5[o][k]] = -1;
                }
                // trivial orbit（只有 1 個要設為 -1）
                seq_b[ orbits1[b1_idx] ] = -1;

                // 檢查 sum(B)：82 - 2*36 = 10
                int sumB = 0;
                for (int i = 0; i < L; i++) sumB += seq_b[i];
                if (sumB != 10) {
                    printf("Warning: B sum != 10 (sumB=%d)\n", sumB);
                }

                // PACF(B)
                for (int u = 0; u < L; u++)
                    pacf_b[u] = periodic_autocorrelation(seq_b, u);

                // 對本批所有 A 檢查
                for (int a = 0; a < batch_count; a++) {
                    checked_pairs++;

                    if (checked_pairs % 1000000ULL == 0ULL) {
                        double percent =
                            (double)A_done * 100.0 / (double)TOTAL_A_COMBOS;
                        printf("Checked %llu pairs so far "
                               "(A combos done: %llu / %llu, %.4f%%).\n",
                               checked_pairs, A_done,
                               TOTAL_A_COMBOS, percent);
                    }

                    if (verify_pcp_with_precomputed(
                            pacfA_batch[a], maxAbsA_batch[a], pacf_b)) {

                        found_pairs++;
                        printf("\n=== Found PCP pair #%llu ===\n", found_pairs);
                        print_sequence(seqA_batch[a], "A");
                        print_sequence(seq_b,          "B");

                        fprintf(fp, "=== Found PCP pair #%llu ===\n", found_pairs);
                        write_sequence(fp, seqA_batch[a], "A");
                        write_sequence(fp, seq_b,          "B");
                        fprintf(fp, "\n");
                        fflush(fp);
                    }
                }

                // 下一個 B 組合
                if (!next_combination(B_comb5, NUM_ORBITS5, B_K5)) {
                    B_has_next = 0;
                }
            } // end while B_has_next
        } // end for each trivial orbit

        double percent_done =
            (double)A_done * 100.0 / (double)TOTAL_A_COMBOS;
        printf("\n--- Finished A batch (A_done = %llu / %llu, %.4f%%) ---\n",
               A_done, TOTAL_A_COMBOS, percent_done);
    }

    printf("\nSearch finished.\n");
    printf("Total A combos processed = %llu / %llu\n",
           A_done, TOTAL_A_COMBOS);
    printf("Total pairs checked      = %llu\n", checked_pairs);
    printf("Total PCP pairs found    = %llu\n", found_pairs);

    fclose(fp);
    return 0;
}