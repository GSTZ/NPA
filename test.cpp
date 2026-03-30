#include <iostream>

using namespace std;
int main() {

    int f, fcmax;
    fcmax = 6;
    f = fcmax % 2 == 0 ? fcmax / 2 : fcmax / 2 + 1;
    cout << f << endl;

    fcmax = 7;
    f = fcmax % 2 == 0 ? fcmax / 2 : fcmax / 2 + 1;
    cout << f << endl;

    return 0;
}


Eigen::MatrixXd ham1(const std::vector<basism>& bas,
                     const SparseMatrix4Dcc& ystrm, // [输入已改为稀疏结构]
                     std::vector<std::map<int, Matrix4D>> buildVValue1,
                     std::vector<double> strength)
{
    auto start = std::chrono::high_resolution_clock::now();

    int nucnum = nucleusm.size();
    int sizebas = bas.size();
    const int t_size = buildVValue1.size();

    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(sizebas, sizebas);

    // 1. 预计算 Interaction Terms (ocal)
    // 这一步生成的是稠密矩阵，保留原样
    std::vector<Matrix4D> ocal_cache(t_size);
    const bool debug_ocal_cache = false;   // 调试开关：是否输出 ocal_cache 信息
    const int debug_ocal_limit = std::max(1, t_size); // 完整明细很大，默认只打印前2个
    const std::string ocal_dump_file = "ocal_cache_dump.txt";
    #pragma omp parallel for schedule(dynamic)
    for (int t = 0; t < t_size; ++t) {
        ocal_cache[t] = calo(buildVValue1[t]);
    }

    auto end_prep = std::chrono::high_resolution_clock::now();
    // std::cout << "Ocal prep done.\n";

    // 2. 主循环
    size_t total_pairs = sizebas * sizebas;
    std::atomic<size_t> completed_pairs(0);
    const size_t update_interval = std::max<size_t>(10, total_pairs / 100);
    auto start_time = std::chrono::steady_clock::now();

    #pragma omp parallel for
    // #pragma omp parallel for schedule(dynamic)
    for (int l = 0; l < sizebas; ++l)
    {
        for (int m = l; m < sizebas; ++m)
        {
            const basism& bas1 = bas[l];
            const basism& bas2 = bas[m];

            // 宇称检查
            if (bas1.parity != bas2.parity) continue;

            // [关键步骤] 适配数据结构
            // 从 SparseMatrix4Dcc 提取 vector<SpMat> 传给 calg_sp
            const SpMat* ptr_l = ystrm[l];
            std::vector<SpMat> ystrm_l(ptr_l, ptr_l + ystrm.cols);

            const SpMat* ptr_m = ystrm[m];
            std::vector<SpMat> ystrm_m(ptr_m, ptr_m + ystrm.cols);

            // [计算 GS] 得到稀疏的 geometric matrix
            SparseMatrix4Dcc gs_sp = calg_sp(bas1, bas2, ystrm_l, ystrm_m);

            // [计算 Hsum] 遍历所有相互作用项
            for (int t = 0; t < t_size; ++t)
            {
                // [核心优化] 替代 multip4dm_sp
                // 直接收缩，不产生中间矩阵
                double hsum = contract_sparse_dense(gs_sp, ocal_cache[t]);

                if (std::abs(hsum) > 1e-15) {
                    // 无锁累加 (不同线程写不同 l,m，安全)
                    result(l, m) += hsum * strength[t];
                }
            }
            result(m, l) = result(l, m);

            // --- 进度打印 (主线程处理，减少锁竞争) ---
            size_t current = ++completed_pairs;
            if (omp_get_thread_num() == 0) {
                if (current % update_interval == 0 || current == total_pairs) {
                    auto now = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
                    std::cout << "\rProcessed: " << current << "/" << total_pairs
                              << " (" << std::fixed << std::setprecision(1) << 100.0 * current / total_pairs << "%)"
                              << " | Time: " << elapsed << "s" << std::flush;
                }
            }
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "\nTotal ham1 time: " << duration.count() << "μs\n";
    if (!outfile.is_open())
    {
        // 如果文件未打开，则尝试打开文件
        outfile.open("basis_output.txt",std::ios::app);
    }
    outfile << "\nTotal ham1 time: " << duration.count() << "μs\n";

    return result;
}