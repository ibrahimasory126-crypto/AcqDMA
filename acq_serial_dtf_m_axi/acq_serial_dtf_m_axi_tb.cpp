#include "acq_serial_dtf_m_axi.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <direct.h>

static void print_csim_cwd() {
    char cwd[2048];
    if (_getcwd(cwd, sizeof(cwd)) != nullptr) {
        std::cout << "[INFO] csim cwd = " << cwd << std::endl;
    } else {
        std::cout << "[WARN] impossible de lire le cwd." << std::endl;
    }
}

static bool parse_int_token(const std::string &token, int &value) {
    const char *s = token.c_str();
    while (*s && std::isspace((unsigned char)*s)) {
        ++s;
    }
    if (*s == '\0') {
        return false;
    }

    char *end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (end == s) {
        return false;
    }

    while (*end && std::isspace((unsigned char)*end)) {
        ++end;
    }
    if (*end != '\0') {
        return false;
    }

    value = (int)v;
    return true;
}

static bool load_csv_int32_one(const std::string &path, std::vector<int> &out) {
    std::ifstream fin(path.c_str());
    if (!fin) {
        return false;
    }

    out.clear();
    std::string line;
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) {
            int v = 0;
            if (parse_int_token(token, v)) {
                out.push_back(v);
            }
        }
    }
    return !out.empty();
}

static bool load_csv_int32_candidates(
    const std::vector<std::string> &candidates,
    const char *label,
    std::vector<int> &out
) {
    for (size_t i = 0; i < candidates.size(); i++) {
        if (load_csv_int32_one(candidates[i], out)) {
            std::cout << "[OK] " << label << " charge depuis: " << candidates[i]
                      << " (n=" << out.size() << ")" << std::endl;
            return true;
        }
    }

    std::cout << "[ERREUR] " << label << " introuvable. chemins testes:" << std::endl;
    for (size_t i = 0; i < candidates.size(); i++) {
        std::cout << "  - " << candidates[i] << std::endl;
    }
    return false;
}

int main() {
    std::cout << "==============================================" << std::endl;
    std::cout << "TB acquisition_serial_m_axi (PRN non-packe)" << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "[INFO] N=" << N
              << " NB_PHASES=" << NB_PHASES
              << " FD_START=" << FD_START
              << " FD_END=" << FD_END << std::endl;

#ifndef __SYNTHESIS__
    std::cout << "[WARN] C-sim utilise la config du header (normal ou CSIM_FAST)." << std::endl;
#endif

    print_csim_cwd();

    std::vector<int> sig_vec;
    std::vector<int> prn_vec;

    std::vector<std::string> sig_paths;
    sig_paths.push_back("signal.csv");
    sig_paths.push_back("../signal.csv");
    sig_paths.push_back("../../signal.csv");
    sig_paths.push_back("../../../signal.csv");

    std::vector<std::string> prn_paths;
    prn_paths.push_back("PRN-25.csv");
    prn_paths.push_back("../PRN-25.csv");
    prn_paths.push_back("../../PRN-25.csv");
    prn_paths.push_back("../../../PRN-25.csv");

    bool ok_sig = load_csv_int32_candidates(sig_paths, "signal.csv", sig_vec);
    bool ok_prn = load_csv_int32_candidates(prn_paths, "PRN-25.csv", prn_vec);

    if (!ok_sig || !ok_prn) {
        std::cerr << "[ECHEC] Entrees indisponibles." << std::endl;
        return 1;
    }

    if ((int)sig_vec.size() < N || (int)prn_vec.size() < N) {
        std::cerr << "[ERREUR] taille insuffisante: signal=" << sig_vec.size()
                  << " prn=" << prn_vec.size() << " N=" << N << std::endl;
        return 1;
    }

    // AXI Full: mem_word_t (normalement ap_int<32>)
    static mem_word_t rx_real[N];
    static mem_word_t prn_in[N];

    const int fd_step = 250;
    const int nb_fd = (FD_END - FD_START) / fd_step + 1;
    const int total_corr = nb_fd * NB_PHASES;

    std::vector<mem_word_t> corr_out_vec(total_corr, 0);
    mem_word_t *corr_out = corr_out_vec.data();

    for (int i = 0; i < N; i++) {
        int s = sig_vec[i];
        if (s == 255) {
            s = -1;
        }

        // PRN brute, meme convention que ton flux precedent
        int p = prn_vec[i];

        rx_real[i] = (s > 0) ? mem_word_t(1) : mem_word_t(-1);
        prn_in[i]  = (mem_word_t)p;
    }

    // Sorties top AXI Full (references int)
    int corr_count     = 0;
    int doppler_out    = 0;
    int codephase_out  = 0;
    int sat_detected   = 0;
    int max_power_out  = 0;
    int mean_power_out = 0;

    int rx_count       = 0;
    int prn_count      = 0;
    int rx_last_seen   = 0;
    int prn_last_seen  = 0;
    int rx_last_pos    = -1;
    int prn_last_pos   = -1;

    acquisition_serial_m_axi(
        rx_real,
        prn_in,
        corr_out,
        corr_count,
        doppler_out,
        codephase_out,
        sat_detected,
        fd_step,
        max_power_out,
        mean_power_out,
        rx_count,
        prn_count,
        rx_last_seen,
        prn_last_seen,
        rx_last_pos,
        prn_last_pos
    );

    bool ok = true;

    if (rx_count != N) {
        std::cerr << "[ERREUR] rx_count=" << rx_count << " attendu=" << N << std::endl;
        ok = false;
    }
    if (prn_count != N) {
        std::cerr << "[ERREUR] prn_count=" << prn_count << " attendu=" << N << std::endl;
        ok = false;
    }
    if (corr_count != total_corr) {
        std::cerr << "[ERREUR] corr_count=" << corr_count << " attendu=" << total_corr << std::endl;
        ok = false;
    }

    int best_idx = 0;
    mem_word_t best_val_sw = corr_out[0];
    for (int i = 1; i < total_corr; i++) {
        if (corr_out[i] > best_val_sw) {
            best_val_sw = corr_out[i];
            best_idx = i;
        }
    }

    int best_fd_idx_sw = best_idx / NB_PHASES;
    int best_tau_sw = best_idx % NB_PHASES;
    int doppler_sw = FD_START + best_fd_idx_sw * fd_step;

    std::cout << "\n[RESULTATS DUT]" << std::endl;
    std::cout << "doppler_out    = " << doppler_out << std::endl;
    std::cout << "codephase_out  = " << codephase_out << std::endl;
    std::cout << "max_power_out  = " << max_power_out << std::endl;
    std::cout << "mean_power_out = " << mean_power_out << std::endl;
    std::cout << "sat_detected   = " << sat_detected << std::endl;
    std::cout << "rx_count       = " << rx_count << std::endl;
    std::cout << "prn_count      = " << prn_count << std::endl;
    std::cout << "corr_count     = " << corr_count << std::endl;
    std::cout << "rx_last_seen   = " << rx_last_seen << std::endl;
    std::cout << "prn_last_seen  = " << prn_last_seen << std::endl;
    std::cout << "rx_last_pos    = " << rx_last_pos << std::endl;
    std::cout << "prn_last_pos   = " << prn_last_pos << std::endl;

    std::cout << "\n[CHECK corr_out]" << std::endl;
    std::cout << "best_idx_sw    = " << best_idx << std::endl;
    std::cout << "best_val_sw    = " << (int)best_val_sw << std::endl;
    std::cout << "doppler_sw     = " << doppler_sw << std::endl;
    std::cout << "phase_sw       = " << best_tau_sw << std::endl;

    if (doppler_out != doppler_sw) {
        std::cerr << "[ERREUR] doppler_out != doppler_sw" << std::endl;
        ok = false;
    }
    if (codephase_out != best_tau_sw) {
        std::cerr << "[ERREUR] codephase_out != phase_sw" << std::endl;
        ok = false;
    }
    if (max_power_out != (int)best_val_sw) {
        std::cerr << "[ERREUR] max_power_out != best_val_sw" << std::endl;
        ok = false;
    }

    if (!ok) {
        std::cerr << "\n[ECHEC] Testbench KO." << std::endl;
        return 1;
    }

    std::cout << "\n[OK] Testbench passe." << std::endl;
    return 0;
}