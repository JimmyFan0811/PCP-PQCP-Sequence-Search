#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
#ifdef _OPENMP

#else
int omp_get_thread_num() { return 0; }
int omp_get_num_threads() { return 1; }
double omp_get_wtime() { return (double)clock() / CLOCKS_PER_SEC; }
#endif

// --- 針對 Project 2 的參數設定 ---
#define L 44               // 目標長度 (可改為 46, 58...)
#define TARGET_MAG 4       // 目標非零值的大小
#define TARGET_COUNT 2     // 目標非零值的數量 (Two nonzero sums)
#define MAX_ITER 5000000   // 單次搜尋最大迭代
#define STAGNATION_LIMIT 40000
#define INIT_TEMP 5.0
#define COOLING_RATE 0.999995

// --- Windows 相容性 ---
#ifdef _MSC_VER
  #include <intrin.h>
  #define POPCOUNT(x) __popcnt64(x)
#else
  #define POPCOUNT(x) __builtin_popcountll(x)
#endif

// --- 全域變數 ---
uint64_t mask;
volatile int solution_found = 0;
// 用來儲存找到的結果以便寫入檔案
uint64_t final_a = 0;
uint64_t final_b = 0;
int success_thread_id = -1;
int success_restart_cnt = 0;
int success_iter = 0;

// --- 高速亂數產生器 ---
typedef struct { uint32_t state; } RngState;
uint32_t xorshift32(RngState *rng) {
    uint32_t x = rng->state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng->state = x;
    return x;
}
int rand_bit(RngState *rng) { return (xorshift32(rng) & 1); }
int rand_range(RngState *rng, int range) { return xorshift32(rng) % range; }
double rand_double(RngState *rng) { return (double)xorshift32(rng) / (double)UINT32_MAX; }

// --- 核心邏輯 ---
uint64_t rotate_left(uint64_t x, int shift) {
    return ((x << shift) | (x >> (L - shift))) & mask;
}

int calculate_cost(uint64_t a, uint64_t b) {
    int penalty = 0;
    int count_target = 0; 

    for (int u = 1; u < L; u++) {
        uint64_t diff_a = a ^ rotate_left(a, u);
        uint64_t diff_b = b ^ rotate_left(b, u);
        
        int rho_a = L - 2 * (int)POPCOUNT(diff_a);
        int rho_b = L - 2 * (int)POPCOUNT(diff_b);
        
        int sum_mag = abs(rho_a + rho_b);
        
        if (sum_mag == 0) {
            continue;
        } else if (sum_mag == TARGET_MAG) {
            count_target++;
        } else {
            penalty += (sum_mag * 10); 
        }
    }

    if (count_target != TARGET_COUNT) {
        penalty += abs(count_target - TARGET_COUNT) * 50;
    }

    return penalty;
}

uint64_t mutate(uint64_t seq, int num_bits, RngState *rng) {
    uint64_t s = seq;
    for(int i=0; i<num_bits; i++) s ^= (1ULL << rand_range(rng, L));
    return s;
}

uint64_t random_seq(RngState *rng) {
    uint64_t s = 0;
    for (int i = 0; i < L; i++) if (rand_bit(rng)) s |= (1ULL << i);
    return s;
}

// 輔助函式：將序列轉為字串 (方便寫入檔案)
void seq_to_string(uint64_t s, char* buffer) {
    for (int i = 0; i < L; i++) {
        buffer[i] = ((s >> i) & 1) ? '1' : '0';
    }
    buffer[L] = '\0';
}

void print_binary(uint64_t s) {
    printf("(");
    for (int i = 0; i < L; i++) printf("%llu", (s >> i) & 1);
    printf(")\n");
}

int main() {
    mask = (1ULL << L) - 1;
    srand((unsigned int)time(NULL));

    // 準備檔名
    char filename[64];
    sprintf(filename, "PQCP_Result_L%d.txt", L);

    printf("=== Project 2 PQCP Search (L=%d) ===\n", L);
    printf("Condition: Exactly %d nonzero sums with magnitude %d.\n", TARGET_COUNT, TARGET_MAG);
    printf("Output will be saved to: %s\n", filename);
    
    #ifdef _OPENMP
    printf("Running on %d threads.\n", omp_get_max_threads());
    #else
    printf("Running single-threaded.\n");
    #endif

    double start_time = omp_get_wtime();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        RngState rng;
        rng.state = (uint32_t)(time(NULL) ^ (tid * 0x5DEECE66) ^ 0x12345678); 
        for(int k=0; k<10; k++) xorshift32(&rng);

        int restart_cnt = 0;
        
        while (!solution_found) {
            uint64_t a = random_seq(&rng);
            uint64_t b = random_seq(&rng);
            int current_cost = calculate_cost(a, b);
            int best_cost_run = current_cost;
            int stagnation = 0;
            double temp = INIT_TEMP;

            for (int iter = 0; iter < MAX_ITER; iter++) {
                if (solution_found) break;

                if (tid == 0 && iter % 10000 == 0) {
                    double elapsed = omp_get_wtime() - start_time;
                    printf("\r[Time: %.1fs] T0 Status -> Cost:%-3d | Best:%-3d | Temp:%.4f", 
                           elapsed, current_cost, best_cost_run, temp);
                    fflush(stdout);
                }

                uint64_t next_a = a;
                uint64_t next_b = b;
                
                if (rand_bit(&rng)) next_a = mutate(a, 1, &rng);
                else                next_b = mutate(b, 1, &rng);

                int next_cost = calculate_cost(next_a, next_b);
                int delta = next_cost - current_cost;

                if (delta <= 0 || (exp(-delta / temp) > rand_double(&rng))) {
                    a = next_a;
                    b = next_b;
                    current_cost = next_cost;
                    
                    if (current_cost < best_cost_run) {
                        best_cost_run = current_cost;
                        stagnation = 0;
                    }
                }

                if (current_cost == 0) {
                    #pragma omp critical
                    {
                        if (!solution_found) {
                            solution_found = 1;
                            // 儲存全域變數，供檔案寫入使用
                            final_a = a;
                            final_b = b;
                            success_thread_id = tid;
                            success_restart_cnt = restart_cnt;
                            success_iter = iter;

                            printf("\n\n>>> SUCCESS by Thread %d <<<\n", tid);
                        }
                    }
                    break;
                }

                stagnation++;
                if (stagnation > STAGNATION_LIMIT) {
                    temp = INIT_TEMP * 0.5; 
                    a = mutate(a, 3, &rng);
                    b = mutate(b, 3, &rng);
                    current_cost = calculate_cost(a, b);
                    stagnation = 0;
                }

                temp *= COOLING_RATE;
                if (temp < 0.001) temp = INIT_TEMP;
            }
            restart_cnt++;
            
            // 如果您希望程式在嘗試一定次數後自動放棄 (例如跑太久都找不到)，
            // 可以取消註解下面這行 (範例：每個執行緒跑完 5 次 restart 就算失敗)
            // if (restart_cnt > 5) break; 
        }
    } // 平行區塊結束

    // --- 寫入檔案邏輯 ---
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Error: Could not open file for writing.\n");
        return 1;
    }

    double total_time = omp_get_wtime() - start_time;

    if (solution_found) {
        char str_a[L + 1], str_b[L + 1];
        seq_to_string(final_a, str_a);
        seq_to_string(final_b, str_b);

        fprintf(fp, "=== PQCP Search Result: SUCCESS ===\n");
        fprintf(fp, "Length (L): %d\n", L);
        fprintf(fp, "Target: Exactly %d nonzero sums with magnitude %d\n", TARGET_COUNT, TARGET_MAG);
        fprintf(fp, "Time Elapsed: %.2f seconds\n", total_time);
        fprintf(fp, "Found by Thread: %d\n", success_thread_id);
        fprintf(fp, "Restart Count: %d, Iteration: %d\n\n", success_restart_cnt, success_iter);
        
        fprintf(fp, "Sequence A: %s\n", str_a);
        fprintf(fp, "Sequence B: %s\n\n", str_b);

        fprintf(fp, "Verification Sums (|rho(a) + rho(b)|):\n[ ");
        printf("Sequence A: %s\n", str_a);
        printf("Sequence B: %s\n", str_b);
        printf("Verification:\n[ ");
        
        for (int u = 1; u < L; u++) {
            int r = abs((L - 2*(int)POPCOUNT(final_a ^ rotate_left(final_a, u))) + 
                        (L - 2*(int)POPCOUNT(final_b ^ rotate_left(final_b, u))));
            fprintf(fp, "%d ", r);
            printf("%d ", r);
        }
        fprintf(fp, "]\n");
        printf("]\n");

        printf("\nResult saved to %s\n", filename);
    } else {
        // 搜尋不到的情況 (通常只有在您設置了最大 restart 次數限制時才會發生)
        // 或者是手動強制停止 (但手動停止通常不會跑到這裡)
        // 為了演示，如果您修改了程式讓它自動停止，就會寫入這些資訊
        
        fprintf(fp, "=== PQCP Search Result: NOT FOUND ===\n");
        fprintf(fp, "Length (L): %d\n", L);
        fprintf(fp, "Condition: %d nonzero sums with magnitude %d\n", TARGET_COUNT, TARGET_MAG);
        fprintf(fp, "Time Elapsed: %.2f seconds\n", total_time);
        fprintf(fp, "Status: Search completed without finding a solution.\n");
        
        printf("\nSearch finished. No solution found. Log saved to %s\n", filename);
    }

    fclose(fp);
    return 0;
}