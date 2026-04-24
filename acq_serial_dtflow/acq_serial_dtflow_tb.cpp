#include "acq_serial_dtflow.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>

using namespace std;

// #define FS 11999000.0
// #define SEUIL_RATIO 2.5

// =========================================================
// Lecture simplifiée et robuste des fichiers CSV
// =========================================================
static vector<int> read_csv(string file) {
    vector<int> data;
    ifstream f(file);
    if (!f.is_open()) {
        cerr << "Erreur: Impossible d'ouvrir " << file << endl;
        exit(1);
    }

    string line;
    while (getline(f, line)) {
        stringstream ss(line);
        string val;
        // Découpe par virgule
        while (getline(ss, val, ',')) {
            if (!val.empty()) {
                try {
                    data.push_back(stoi(val));
                } catch (...) {
                    // Ignore les entêtes ou caractères non numériques
                }
            }
        }
    }
    f.close();
    return data;
}

static vector<int> read_csv_with_fallback(const vector<string> &paths) {
    for (size_t i = 0; i < paths.size(); i++) {
        ifstream test(paths[i].c_str());
        if (test.is_open()) {
            test.close();
            cout << "[INFO] CSV trouvé : " << paths[i] << endl;
            return read_csv(paths[i]);
        }
    }

    cerr << "[ERREUR] Aucun chemin CSV valide trouvé." << endl;
    for (size_t i = 0; i < paths.size(); i++) {
        cerr << "         - " << paths[i] << endl;
    }
    exit(1);
}

// =========================================================
// MAIN TESTBENCH - CSIM pour version DATAFLOW
// =========================================================
int main() {
    cout << "==================================================" << endl;
    cout << "   TESTBENCH DATAFLOW - CSIM" << endl;
    cout << "==================================================" << endl;
    cout << "N         : " << N << endl;
    cout << "NB_PHASES : " << NB_PHASES << endl;
    cout << "FD_START  : " << FD_START << " Hz" << endl;
    cout << "FD_END    : " << FD_END << " Hz" << endl;
    cout << "FS        : " << FS << " Hz" << endl;
    cout << "==================================================" << endl << endl;

    // =========== Lecture des fichiers CSV ===========
    cout << "[INFO] Lecture des fichiers CSV..." << endl;
    vector<int> signal = read_csv_with_fallback({
        "../signal.csv",
        "../../STSA/signal.csv",
        "C:/Users/LABO01/Desktop/ibra_vscode/STSA/signal.csv"
    });
    vector<int> prn = read_csv_with_fallback({
        "../PRN-25.csv",
        "../../STSA/PRN-25.csv",
        "C:/Users/LABO01/Desktop/ibra_vscode/STSA/PRN-25.csv"
    });

    cout << "[INFO] Signal : " << signal.size() << " samples" << endl;
    cout << "[INFO] PRN    : " << prn.size() << " chips" << endl;

    if (signal.size() < N || prn.size() < N) {
        cerr << "[ERREUR] Taille CSV insuffisante !" << endl;
        cerr << "         Signal: " << signal.size() << " < N=" << N << endl;
        cerr << "         PRN:    " << prn.size() << " < N=" << N << endl;
        return 1;
    }

    // =========== Calcul du nombre de Dopplers ===========
    int fd_step = 250;  // Pas par défaut
    int nb_fd = (FD_END - FD_START) / fd_step + 1;
    int total_outputs = nb_fd * NB_PHASES;

    cout << "[INFO] Pas Doppler (fd_step) : " << fd_step << " Hz" << endl;
    cout << "[INFO] Nombre de Dopplers : " << nb_fd << endl;
    cout << "[INFO] Nombre total de résultats attendus : " << total_outputs << endl << endl;

    // =========== Création des streams ===========
    hls::stream<axis_t> rx_stream;
    hls::stream<axis_t> prn_stream;
    hls::stream<axis_t> corr_stream;

    // =========== Injection des données dans les streams ===========
    cout << "[INFO] Injection des signaux..." << endl;
    for (int i = 0; i < N; i++) {
        axis_t a, b;
        
        // Injection signal
        a.data = signal[i];
        a.last = (i == N - 1) ? 1 : 0;
        a.keep = -1;
        a.strb = -1;
        rx_stream.write(a);
        
        // Injection PRN
        b.data = prn[i];
        b.last = (i == N - 1) ? 1 : 0;
        b.keep = -1;
        b.strb = -1;
        prn_stream.write(b);
    }
    cout << "[INFO] Injection complétée." << endl << endl;

    // =========== Appel de la fonction HLS ===========
    cout << "[INFO] Lancement de acquisition_serial()..." << endl;
    int doppler, codephase, sat;
    int max_power_dbg = 0;    
    int mean_power_dbg = 0;   
    int rx_count_dbg = 0;
    int prn_count_dbg = 0;
    int rx_last_seen_dbg = 0;
    int prn_last_seen_dbg = 0;
    int rx_last_pos_dbg = -1;
    int prn_last_pos_dbg = -1;

    auto start = chrono::high_resolution_clock::now();

    acquisition_serial(
        rx_stream,
        prn_stream,
        corr_stream,
        doppler,
        codephase,
        sat,
        fd_step,
        max_power_dbg,      
        mean_power_dbg,     
        rx_count_dbg,
        prn_count_dbg,
        rx_last_seen_dbg,
        prn_last_seen_dbg,
        rx_last_pos_dbg,
        prn_last_pos_dbg
    );
    
    auto stop = chrono::high_resolution_clock::now();
    double elapsed_ms = chrono::duration<double, milli>(stop - start).count();

    cout << "[INFO] Exécution terminée en " << fixed << setprecision(2) << elapsed_ms << " ms" << endl << endl;

    // =========== Lecture des résultats ===========
    cout << "[INFO] Lecture des résultats de corrélation..." << endl;
    
    power_t max_power = 0;
    int max_fd_idx = 0;
    int max_tau = 0;
    power_t sum_power = 0;
    int count = 0;
    int last_count = 0;
    int last_index = -1;

    int fd_idx_current = 0;
    int tau_current = 0;

    // Buffer pour affichage structuré (optionnel)
    bool verbose = false; // Mettez à true pour affichage détaillé

    while (!corr_stream.empty() && count < total_outputs) {
        axis_t out = corr_stream.read();
        
        int current_power_int = (int)out.data;
        power_t current_power = (power_t)current_power_int;
        
        sum_power += current_power;
        
        // Déterminer fd_idx et tau depuis count
        fd_idx_current = count / NB_PHASES;
        tau_current = count % NB_PHASES;

        // Trouver le pic
        if (current_power > max_power) {
            max_power = current_power;
            max_fd_idx = fd_idx_current;
            max_tau = tau_current;
        }

        if (out.last) {
            last_count++;
            last_index = count;
        }

        count++;

        // Affichage optionnel (très verbeux avec NB_PHASES=1023 et nb_fd=21)
        if (verbose && (count <= 10 || out.last)) {
            cout << "  [" << count << "] fd_idx=" << fd_idx_current 
                 << " tau=" << tau_current 
                 << " power=" << (int)current_power 
                 << " last=" << (int)out.last << endl;
        }
    }

    cout << "[INFO] Nombre de résultats lus : " << count << endl;
    cout << "[INFO] last_count : " << last_count << endl;
    cout << "[INFO] last_index : " << last_index << endl;

    if (count != total_outputs) {
        cerr << "[AVERTISSEMENT] Nombre de résultats inattendus !" << endl;
        cerr << "                Attendu: " << total_outputs << ", reçu: " << count << endl;
    } else {
        cout << "[OK] Tous les résultats reçus." << endl;
    }

    // =========== Calcul des statistiques ===========
    power_t mean_power = 0;
    if (count > 0) {
        mean_power = sum_power / count;
    }

    // =========== Conversion des indices Doppler/Phase ===========
    int doppler_estimated = FD_START + max_fd_idx * fd_step;
    int codephase_estimated = max_tau;

    // =========== Affichage des résultats ===========
    cout << endl;
    cout << "=" << string(48, '=') << endl;
    cout << "                    RESULTATS CSIM" << endl;
    cout << "=" << string(48, '=') << endl;
    cout << endl;
    
    cout << "[SORTIE REGISTRES (AXI-Lite)]" << endl;
    cout << "  Doppler détecté   : " << setw(6) << doppler << " Hz (Attendu: ~-1750)" << endl;
    cout << "  Phase de code     : " << setw(6) << codephase << " (Attendu: ~150)" << endl;
    cout << "  Satellite         : " << (sat ? "DETECTE" : "NON DETECTE") << endl;

    cout << "  Max power (HW)    : " << max_power_dbg << endl;
    cout << "  Mean power (HW)   : " << mean_power_dbg << endl;
    cout << "  Ratio HW max/mean : "
     << (mean_power_dbg > 0 ? (double)max_power_dbg / (double)mean_power_dbg : 0.0)
     << endl;
    
    cout << endl;
    cout << "[RESULTATS CALCULES]" << endl;
    cout << "  Doppler (pic max) : " << setw(6) << doppler_estimated << " Hz" << endl;
    cout << "  Phase (pic max)   : " << setw(6) << codephase_estimated << endl;
    cout << "  Puissance max     : " << (int)max_power << endl;
    cout << "  Puissance moyenne : " << (int)mean_power << endl;
    cout << "  Ratio max/mean    : " << fixed << setprecision(2) 
         << (mean_power > 0 ? (double)max_power / (double)mean_power : 0.0) << endl;
        power_t threshold = mean_power * (power_t)(1.0 + SEUIL_K);
        cout << "  Threshold CFAR (x" << (1.0 + SEUIL_K) << "): " << (int)threshold << endl;
    
    cout << endl;
    cout << "[STATISTIQUES]" << endl;
    cout << "  Points de corrélation lus : " << count << endl;
    cout << "  Temps CSIM total          : " << fixed << setprecision(2) 
         << elapsed_ms << " ms" << endl;
    
    cout << endl;
    cout << "=" << string(48, '=') << endl;
    
    // =========== Vérifications ===========
    if (count == total_outputs) {
        cout << "[PASS] Nombre de résultats correct" << endl;
    } else {
        cout << "[FAIL] Erreur sur le nombre de résultats" << endl;
    }

    if (last_count == 1 && last_index == count - 1) {
        cout << "[PASS] Protocole AXIS last correct" << endl;
    } else {
        cout << "[FAIL] Protocole AXIS last incorrect" << endl;
    }

    if (sat) {
        cout << "[PASS] Satellite détecté sur le seuil" << endl;
    } else {
        cout << "[FAIL] Satellite non détecté (à vérifier)" << endl;
    }

    cout << "=" << string(48, '=') << endl << endl;

    return 0;
}







