#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

// 顯式定義 M_PI
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// =================================================================
// 參數設定
// =================================================================
#define INPUT_FILENAME "Unique_PCPs_82_Processed.txt"
#define OUTPUT_FILENAME "papr_combined_output_len82_Processed.txt" // 檔名稍微修改以區分

#define MAX_SEQ_LEN 256
#define LINE_BUF_SIZE 2048 
#define MAX_SEQ_STR_LEN 1024

// 模擬點數 (Oversampling)
#define SIMULATION_POINTS 2000 

// 儲存序列結構
typedef struct {
    double data[MAX_SEQ_LEN];
    int length;
    char original_str[MAX_SEQ_STR_LEN];
} Sequence;

// =================================================================
// 核心計算函式 (Linear PAPR)
// =================================================================
double calculate_papr_linear(Sequence *seq) {
    int N = seq->length;
    if (N == 0) return 999.0;

    double max_inst_power = 0.0;
    double total_inst_power = 0.0;

    for (int i = 0; i < SIMULATION_POINTS; i++) {
        double t = (double)i / (double)SIMULATION_POINTS;
        double complex S = 0.0 + 0.0 * I;

        for (int n = 0; n < N; n++) {
            double angle = 2.0 * M_PI * t * n;
            S += seq->data[n] * (cos(angle) + I * sin(angle));
        }
        
        S = S * (1.0 / sqrt((double)N));

        double p_inst = creal(S) * creal(S) + cimag(S) * cimag(S);

        if (p_inst > max_inst_power) {
            max_inst_power = p_inst;
        }
        total_inst_power += p_inst;
    }

    double p_avg = total_inst_power / SIMULATION_POINTS;
    if (p_avg == 0) return 0.0;
    
    return max_inst_power / p_avg;
}

// =================================================================
// 解析函式
// =================================================================
int parse_sequence(const char *line, const char *prefix, Sequence *out_seq) {
    char *start_ptr = strstr(line, prefix);
    if (!start_ptr) return 0;

    char *brace_start = strchr(start_ptr, '{'); 
    char *brace_end = strchr(brace_start, '}');   
    
    if (!brace_start || !brace_end) return 0;

    int len_str = brace_end - brace_start + 1;
    if (len_str >= MAX_SEQ_STR_LEN) len_str = MAX_SEQ_STR_LEN - 1;
    
    strncpy(out_seq->original_str, brace_start, len_str);
    out_seq->original_str[len_str] = '\0'; 

    out_seq->length = 0;
    char *current_char = brace_start + 1; 

    while (*current_char != '\0' && *current_char != '}') {
        if (*current_char == '+') {
            out_seq->data[out_seq->length++] = 1.0;
        } else if (*current_char == '-') {
            out_seq->data[out_seq->length++] = -1.0;
        }
        current_char++;
    }

    return 1;
}

// =================================================================
// 主程式
// =================================================================
int main() {
    FILE *fin = fopen(INPUT_FILENAME, "r");
    FILE *fout = fopen(OUTPUT_FILENAME, "w");

    if (!fin) {
        printf("Error: Cannot open input file '%s'\n", INPUT_FILENAME);
        return 1;
    }
    if (!fout) {
        printf("Error: Cannot open output file '%s'\n", OUTPUT_FILENAME);
        fclose(fin);
        return 1;
    }

    char line[LINE_BUF_SIZE];
    Sequence seqA, seqB;
    int foundA = 0, foundB = 0;
    int pair_idx = 0;
    
    // 全局最佳紀錄 (這次我們要找 Min-Max)
    // 意思是：找出 "該 Pair 較大的那個 PAPR" 最小的那一組
    double global_min_max_papr = 1e9;
    int best_pair_idx = -1;
    
    // 記錄造成該 Pair 成為 Max 的那個序列 (瓶頸序列)
    char bottleneck_seq_str[MAX_SEQ_STR_LEN] = "";
    char bottleneck_seq_type = ' '; 

    printf("Processing %s (Min-Max PAPR Strategy)...\n", INPUT_FILENAME);
    fprintf(fout, "PAPR Analysis Results (Min-Max Strategy)\n");
    fprintf(fout, "Strategy: Minimize the Max(PAPR_A, PAPR_B) across all pairs.\n");
    fprintf(fout, "Input: %s\n", INPUT_FILENAME);
    fprintf(fout, "---------------------------------------------------------------------------------------------------\n");

    while (fgets(line, sizeof(line), fin)) {
        char *pair_ptr = strstr(line, "Pair #");
        if (pair_ptr != NULL) {
            foundA = 0;
            foundB = 0;
            sscanf(pair_ptr, "Pair #%d", &pair_idx);
        }

        if (parse_sequence(line, "A =", &seqA)) foundA = 1;
        if (parse_sequence(line, "B =", &seqB)) foundB = 1;

        if (foundA && foundB) {
            double paprA_lin = calculate_papr_linear(&seqA);
            double paprB_lin = calculate_papr_linear(&seqB);
            
            // ========================================================
            // 修改重點：找出此 Pair 中的 "最大值" (Worst case in this pair)
            // ========================================================
            double pair_max_papr_lin = (paprA_lin > paprB_lin) ? paprA_lin : paprB_lin;
            double pair_max_papr_db = 10.0 * log10(pair_max_papr_lin);

            // ========================================================
            // 修改重點：在所有 Pair 的 "最大值" 中，找出 "最小值"
            // (我們希望最差的情況越低越好)
            // ========================================================
            if (pair_max_papr_lin < global_min_max_papr) {
                global_min_max_papr = pair_max_papr_lin;
                best_pair_idx = pair_idx;
                
                // 記錄是誰決定了這個 Pair 的 Max 值
                if (paprA_lin > paprB_lin) {
                    strcpy(bottleneck_seq_str, seqA.original_str);
                    bottleneck_seq_type = 'A';
                } else {
                    strcpy(bottleneck_seq_str, seqB.original_str);
                    bottleneck_seq_type = 'B';
                }
            }
            
            // 寫入檔案 (標示 Pair Max)
            fprintf(fout, "\n=== Pair #%d (Max PAPR: %.4f dB / %.4f linear) ===\n", 
                    pair_idx, pair_max_papr_db, pair_max_papr_lin);
            fprintf(fout, "A (PAPR=%.4f): %s\n", paprA_lin, seqA.original_str);
            fprintf(fout, "B (PAPR=%.4f): %s\n", paprB_lin, seqB.original_str);
            
            if (pair_idx % 10 == 0) {
                printf("Processed Pair #%d (Current Min-Max: %.4f dB)\n", 
                       pair_idx, 10.0 * log10(global_min_max_papr));
            }

            foundA = 0;
            foundB = 0;
        }
    }

    fprintf(fout, "\n===================================================\n");
    fprintf(fout, "BEST RESULT FOUND (Min-Max Strategy)\n");
    fprintf(fout, "===================================================\n");
    fprintf(fout, "Pair Index: #%d\n", best_pair_idx);
    fprintf(fout, "Lowest 'Pair-Max' PAPR (Linear): %.6f\n", global_min_max_papr);
    fprintf(fout, "Lowest 'Pair-Max' PAPR (dB):     %.6f dB\n", 10.0 * log10(global_min_max_papr));
    fprintf(fout, "Bottleneck Sequence (%c): %s\n", bottleneck_seq_type, bottleneck_seq_str);

    printf("\nDone! Best Min-Max PAPR found: %.6f dB (Pair #%d).\n", 
           10.0 * log10(global_min_max_papr), best_pair_idx);
    printf("Results saved to: %s\n", OUTPUT_FILENAME);

    fclose(fin);
    fclose(fout);

    return 0;
}