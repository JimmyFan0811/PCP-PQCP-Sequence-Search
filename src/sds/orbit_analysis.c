#include <stdio.h>
#include <stdlib.h>
#include <string.h> // For memset
#include <math.h>   // For sqrt

// 設定一個合理的陣列大小上限
#define MAX_V 1000 
#define MAX_DIVISORS 1000

// ==================================================
// 輔助函數
// ==================================================

// 計算最大公因數 (GCD)
int gcd(int a, int b) {
    while (b) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// 計算 (base^exp) % mod (模冪運算)
long long power(long long base, int exp, int mod) {
    long long res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return res;
}

// ==================================================
// 步驟 0: 找出可能的 (r, s) 配對
// (來自 05_PCP.pdf, Theorem 1) [cite: 259, 260, 262]
// ==================================================

// 檢查一個數是否為平方數
int is_perfect_square(int n) {
    if (n < 0) return 0;
    int root = (int)(sqrt(n) + 0.5);
    return root * root == n;
}

// 找出 (r, s) 候選配對
void find_r_s_pairs(int v) {
    printf("==================================================\n");
    printf("步驟 0: 尋找 (r, s) 候選配對 (基於 Theorem 1)\n");
    printf("       (L-2r)^2 + (L-2s)^2 = 2L (L=%d)\n", v);
    printf("       尋找 (A')^2 + (B')^2 = %d\n\n", 2 * v);

    int target = 2 * v;
    int found = 0;

    // 迭代 A' 從 0 到 sqrt(target)
    for (int a_prime = 0; a_prime * a_prime <= target; a_prime++) {
        int b_prime_sq = target - (a_prime * a_prime);
        
        if (is_perfect_square(b_prime_sq)) {
            int b_prime = (int)(sqrt(b_prime_sq) + 0.5);
            
            // 我們找到了一對 (A', B') = (a_prime, b_prime)
            // 現在反解 r 和 s
            
            // A' = v - 2r  => 2r = v - a_prime
            if ((v - a_prime) % 2 == 0) {
                int r1 = (v - a_prime) / 2;
                // B' = v - 2s  => 2s = v - b_prime
                if ((v - b_prime) % 2 == 0) {
                    int s1 = (v - b_prime) / 2;
                    printf("  找到一組解: A' = %d, B' = %d  => (r, s) = (%d, %d)\n", a_prime, b_prime, r1, s1);
                    found = 1;
                }
                // B' = -(v - 2s) => 2s = v + b_prime
                if ((v + b_prime) % 2 == 0) {
                    int s2 = (v + b_prime) / 2;
                    printf("  找到一組解: A' = %d, B' = -%d => (r, s) = (%d, %d)\n", a_prime, b_prime, r1, s2);
                    found = 1;
                }
            }

            // A' = -(v - 2r) => 2r = v + a_prime
            if (a_prime != 0 && (v + a_prime) % 2 == 0) {
                int r2 = (v + a_prime) / 2;
                // B' = v - 2s  => 2s = v - b_prime
                if ((v - b_prime) % 2 == 0) {
                    int s1 = (v - b_prime) / 2;
                    printf("  找到一組解: A' = -%d, B' = %d  => (r, s) = (%d, %d)\n", a_prime, b_prime, r2, s1);
                    found = 1;
                }
                // B' = -(v - 2s) => 2s = v + b_prime
                if ((v + b_prime) % 2 == 0) {
                    int s2 = (v + b_prime) / 2;
                    printf("  找到一組解: A' = -%d, B' = -%d => (r, s) = (%d, %d)\n", a_prime, b_prime, r2, s2);
                    found = 1;
                }
            }
        }
    }
    
    if (!found) {
        printf("  未找到滿足 Theorem 1 的 (r, s) 整數解。\n");
    }
}


// ==================================================
// 步驟 1: 計算 phi(v) (歐拉函數)
// 筆記中 size=36 的計算
// ==================================================
int phi(int v) {
    int count = 0;
    for (int i = 1; i < v; i++) {
        if (gcd(v, i) == 1) {
            count++;
        }
    }
    return count;
}

// ==================================================
// 步驟 2 & 3: 軌道分析
// ==================================================

// 輔助函數：計算一個元素 g 在 Z_v^* 中的「階」 (Order)
int get_order(int g, int v) {
    if (g <= 0 || gcd(g, v) != 1) {
        return -1; // 不在 Z_v^* 中
    }
    
    long long p = g;
    for (int k = 1; k < v; k++) {
        if (p == 1) {
            return k; // 階是 k
        }
        p = (p * g) % v;
    }
    return -1; 
}

// 找出 n 的所有因數
void find_divisors(int n, int divisors[], int* count) {
    *count = 0;
    for (int i = 1; i <= n; i++) {
        if (n % i == 0) {
            divisors[(*count)++] = i;
        }
    }
}

// 尋找並印出軌道
void find_orbits(int v, int H_elements[], int H_size) {
    printf("  Orbits of Z_%d under H:\n", v);
    
    int visited[MAX_V];
    memset(visited, 0, sizeof(visited)); // 初始化為 0
    
    int total_trivial = 0;
    int total_nontrivial = 0;

    // 遍歷 Z_v 中的每一個元素 j
    //
    for (int j = 0; j < v; j++) {
        if (!visited[j]) {
            int orbit[MAX_V];
            int orbit_size = 0;
            
            // 計算 H * j
            for (int k = 0; k < H_size; k++) {
                int h = H_elements[k];
                int orbit_element = (1LL * h * j) % v;
                
                int found = 0;
                for (int m = 0; m < orbit_size; m++) {
                    if (orbit[m] == orbit_element) {
                        found = 1;
                        break;
                    }
                }
                
                if (!found) {
                    orbit[orbit_size++] = orbit_element;
                    visited[orbit_element] = 1;
                }
            }
            
            // 印出軌道
            printf("    { ");
            for (int m = 0; m < orbit_size; m++) {
                printf("%d ", orbit[m]);
            }
            printf("}");
            
            // 判斷 Trivial 還是 Nontrivial
            //
            if (orbit_size == 1) {
                printf(" (Trivial)\n");
                total_trivial++;
            } else {
                printf(" (Nontrivial)\n");
                total_nontrivial++;
            }
        }
    }
    printf("  Total Trivial Orbits: %d, Total Nontrivial Orbits: %d\n", total_trivial, total_nontrivial);
}

// ==================================================
// 主程式
// ==================================================

int main() {
    int v;
    printf("請輸入序列長度 (v): ");
    scanf("%d", &v);

    if (v <= 0 || v > MAX_V) {
        printf("v 必須大於 0 且小於 %d\n", MAX_V);
        return 1;
    }

    // 步驟 0: 找出 (r, s) 配對
    find_r_s_pairs(v);

    printf("\n==================================================\n");

    // 步驟 1: 計算 phi(v)
    int phiv = phi(v);
    printf("步驟 1: 乘法群 Z_%d^* 的大小 (phi(%d)) = %d\n", v, v, phiv);

    // 步驟 2: 找出 phi(v) 的所有因數 (可能的 order)
    int divisors[MAX_DIVISORS];
    int divisor_count = 0;
    find_divisors(phiv, divisors, &divisor_count);
    
    printf("步驟 2: Z_%d^* 中可能的子群階 (Order) (phi(%d) 的因數) 有 %d 個:\n  { ", v, phiv, divisor_count);
    for (int i = 0; i < divisor_count; i++) {
        printf("%d ", divisors[i]);
    }
    printf("}\n");

    // 步驟 3: 對每一個 order，找出 H 及其軌道
    for (int i = 0; i < divisor_count; i++) {
        int d = divisors[i];
        printf("\n--------------------------------------------------\n");
        printf("步驟 3: 分析 Order d = %d\n", d);

        // 3a: 尋找一個階為 d 的生成元 g
        int g = -1;
        for (int j = 1; j < v; j++) {
            if (get_order(j, v) == d) {
                g = j;
                break; // 找到一個就夠了
            }
        }

        if (g == -1) {
            printf("  未能在 Z_%d^* 中找到階為 %d 的生成元。\n", v, d);
            continue;
        }

        // 3b: 構造 H 並列出元素
        printf("  3b: 找到一個生成元 g = %d\n", g);
        printf("      對應的循環子群 H = { ");
        
        int H_elements[MAX_DIVISORS];
        for (int k = 0; k < d; k++) {
            H_elements[k] = (int)power(g, k, v);
            printf("%d ", H_elements[k]);
        }
        printf("}\n\n");

        // 3c: 找出 Trivial 和 Nontrivial 軌道
        printf("  3c: 找出 H 作用在 Z_%d 上的軌道:\n", v);
        find_orbits(v, H_elements, d);
    }

    return 0;
}