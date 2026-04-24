#include "acq_serial_dtflow.h"
#include "dds_lut_rom.h"

// Nombre de phases de code traitées en parallèle par tile
static const int TAU_TILE = 8;

// ============================================================
// ÉTAPE 1: LECTURE DES ENTRÉES AXI-STREAM
// ============================================================

/**
 * Lit N échantillons depuis les streams AXI et les pousse dans des FIFOs internes.
 * Convertit le PRN en signe binaire (1 si >0, 0 sinon).
 */
static void read_inputs_to_stream(
    hls::stream<axis_t> &rx_real,
    hls::stream<axis_t> &prn_in,
    hls::stream<data_t> &sig_s,
    hls::stream<ap_uint<1> > &prn_s,
    int &rx_last_seen,
    int &prn_last_seen,
    int &rx_last_pos,
    int &prn_last_pos
) {
#pragma HLS INLINE off
#if !ACQ_ENABLE_TLAST_DEBUG
    (void)rx_last_seen;
    (void)prn_last_seen;
    (void)rx_last_pos;
    (void)prn_last_pos;
#endif
    READ_INPUTS: for (ap_uint<14> i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        axis_t rx  = rx_real.read();
        axis_t prn = prn_in.read();

        sig_s.write((data_t)rx.data);
        // Binarisation PRN: +1 → 1, -1 → 0
        prn_s.write((prn.data > 0) ? (ap_uint<1>)1 : (ap_uint<1>)0);

#if ACQ_ENABLE_TLAST_DEBUG
        // Capture position du signal TLAST pour debug
        if (rx.last) {
            rx_last_seen = 1;
            rx_last_pos  = (int)i;
        }
        if (prn.last) {
            prn_last_seen = 1;
            prn_last_pos  = (int)i;
        }
#endif
    }
}

/**
 * Transfert des FIFOs vers les buffers mémoire locaux (BRAM).
 */
static void store_stream_to_mem(
    hls::stream<data_t> &sig_s,
    hls::stream<ap_uint<1> > &prn_s,
    data_t signal_buf[N],
    ap_uint<1> prn_sign[2 * N]
) {
#pragma HLS INLINE off

    STORE_INPUTS: for (ap_uint<14> i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        signal_buf[i] = sig_s.read();
        prn_sign[i]   = prn_s.read();
    }
}

/**
 * Région dataflow: lecture et stockage en parallèle via FIFOs.
 */
static void capture_inputs_to_mem(
    hls::stream<axis_t> &rx_real,
    hls::stream<axis_t> &prn_in,
    data_t signal_buf[N],
    ap_uint<1> prn_sign[2 * N],
    int &rx_last_seen,
    int &prn_last_seen,
    int &rx_last_pos,
    int &prn_last_pos
) {
#pragma HLS INLINE off

    // FIFOs de découplage pour dataflow
    hls::stream<data_t> sig_s("sig_s");
    hls::stream<ap_uint<1> > prn_s("prn_s");
#pragma HLS STREAM variable=sig_s depth=64
#pragma HLS STREAM variable=prn_s depth=64

#pragma HLS DATAFLOW
    read_inputs_to_stream(rx_real, prn_in, sig_s, prn_s, rx_last_seen, prn_last_seen, rx_last_pos, prn_last_pos);
    store_stream_to_mem(sig_s, prn_s, signal_buf, prn_sign);
}

// ============================================================
// UTILITAIRE DDS
// ============================================================

/**
 * Convertit une fréquence (Hz) en incrément de phase pour le NCO.
 * Formule: phase_inc = freq * 2^32 / Fs (avec arrondi)
 */
static inline ap_int<32> hz_to_phase_inc(int freq_hz) {
#pragma HLS INLINE
    const long long SCALE = (1LL << DDS_PHASE_BITS);
    long long num = (long long)freq_hz * SCALE;
    // Arrondi au plus proche
    if (num >= 0) num += (FS_INT / 2);
    else          num -= (FS_INT / 2);
    return (ap_int<32>)(num / (long long)FS_INT);
}

// ============================================================
// ÉTAPE 2: PRÉPARATION DES DONNÉES
// ============================================================

/**
 * Charge les entrées et prépare le buffer PRN étendu (2×N pour décalages circulaires).
 */
static void load_inputs_once(
    hls::stream<axis_t> &rx_real,
    hls::stream<axis_t> &prn_in,
    data_t signal_buf[N],
    ap_uint<1> prn_sign[2 * N],
    int &rx_count,
    int &prn_count,
    int &rx_last_seen,
    int &prn_last_seen,
    int &rx_last_pos,
    int &prn_last_pos
) {
#pragma HLS INLINE off

#if !ACQ_ENABLE_TLAST_DEBUG
    (void)rx_last_seen;
    (void)prn_last_seen;
    (void)rx_last_pos;
    (void)prn_last_pos;
#endif

    rx_count = 0;
    prn_count = 0;

    capture_inputs_to_mem(
        rx_real, prn_in, signal_buf, prn_sign,
        rx_last_seen, prn_last_seen, rx_last_pos, prn_last_pos
    );

    rx_count = N;
    prn_count = N;

    // Duplication du PRN pour permettre décalages circulaires sans modulo
    DUP_PRN: for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        prn_sign[i + N] = prn_sign[i];
    }
}

/**
 * Précalcule l'index de départ dans le buffer PRN pour chaque phase tau.
 * tau_start[tau] = tau * N / NB_PHASES
 */
static void build_tau_start_table(int tau_start_tbl[NB_PHASES]) {
#pragma HLS INLINE off

    BUILD_TAU_START: for (int tau = 0; tau < NB_PHASES; tau++) {
#pragma HLS PIPELINE II=1
        tau_start_tbl[tau] = (int)(((long long)tau * (long long)N) / NB_PHASES);
    }
}

/**
 * Réplique le buffer PRN en TAU_TILE copies pour accès parallèle.
 * Chaque "bank" permet de lire une phase tau différente simultanément.
 */
static void build_prn_banks(
    const ap_uint<1> prn_sign[2 * N],
    ap_uint<1> prn_banks[TAU_TILE][2 * N]
) {
#pragma HLS INLINE off

    BUILD_PRN_BANKS: for (int i = 0; i < 2 * N; i++) {
#pragma HLS PIPELINE II=1
        COPY_BANKS: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS UNROLL
            prn_banks[k][i] = prn_sign[i];
        }
    }
}

// ============================================================
// ÉTAPE 3: MÉLANGEUR NCO (DOWN-CONVERSION DOPPLER)
// ============================================================

/**
 * Mélange le signal avec un oscillateur local à fréquence (Fc + fd).
 * Produit les composantes I et Q (baseband complexe).
 * Utilise une LUT sin/cos pour le NCO.
 */
static void doppler_mixer_to_mem(
    const data_t signal_buf[N],
    osc_t mixI_loc[N],
    osc_t mixQ_loc[N],
    int fd_hz
) {
#pragma HLS INLINE off

    // Incrément de phase négatif pour down-conversion
    const ap_int<32> phase_inc = hz_to_phase_inc(-(FREQUENCE_CENTRALE_HZ + fd_hz));
    ap_uint<32> phase_acc = 0;

    MIX_ALL: for (int n = 0; n < N; n++) {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=mixI_loc inter false
#pragma HLS DEPENDENCE variable=mixQ_loc inter false

        // Index LUT = bits de poids fort de l'accumulateur
        ap_uint<DDS_LUT_BITS> lut_idx =
            phase_acc.range(DDS_PHASE_BITS - 1, DDS_PHASE_BITS - DDS_LUT_BITS);

        osc_t c = DDS_COS_LUT[(int)lut_idx];  // cos(phase)
        osc_t s = DDS_SIN_LUT[(int)lut_idx];  // sin(phase)
        osc_t x = (osc_t)signal_buf[n];

        // Mélange: I = x*cos, Q = -x*sin (convention down-conversion)
        mixI_loc[n] = x * c;
        mixQ_loc[n] = -(x * s);
        
        phase_acc = (ap_uint<32>)(phase_acc + (ap_uint<32>)phase_inc);
    }
}

// ============================================================
// ÉTAPE 4: CORRÉLATION (TAU_TILE PHASES EN PARALLÈLE)
// ============================================================

/**
 * Calcule la corrélation pour TAU_TILE phases de code consécutives.
 * Utilise un accumulateur pair/impair pour éviter les dépendances de boucle.
 */
static void process_tau_tile_v9(
    const osc_t mixI_loc[N],
    const osc_t mixQ_loc[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<pow_pkt_t> &pow_out,
    int fd_idx,
    int tau_base
) {
#pragma HLS INLINE off

    // Accumulateurs I/Q pour chaque tau du tile, séparés pair/impair
    acc_t accI[TAU_TILE][2];
    acc_t accQ[TAU_TILE][2];
#pragma HLS ARRAY_PARTITION variable=accI complete dim=1
#pragma HLS ARRAY_PARTITION variable=accI complete dim=2
#pragma HLS ARRAY_PARTITION variable=accQ complete dim=1
#pragma HLS ARRAY_PARTITION variable=accQ complete dim=2

    // Initialisation des accumulateurs
    INIT_TILE: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS UNROLL
        accI[k][0] = 0; accI[k][1] = 0;
        accQ[k][0] = 0; accQ[k][1] = 0;
    }

    // Boucle principale: accumulation sur N échantillons
    SAMPLE_LOOP: for (int n = 0; n < N; n++) {
#pragma HLS PIPELINE II=1
        osc_t base_i = mixI_loc[n];
        osc_t base_q = mixQ_loc[n];
        bool odd = (n & 1) != 0;  // Alternance pair/impair

        // Traitement parallèle de TAU_TILE phases
        TAU_TILE_LOOP: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS UNROLL
            int tau = tau_base + k;
            if (tau < NB_PHASES) {
                osc_t i = base_i;
                osc_t q = base_q;
                
                // Lecture du signe PRN décalé de tau
                ap_uint<1> s = prn_banks[k][tau_start_tbl[tau] + n];

                // Multiplication par PRN: si s=0 (→-1), inverser le signe
                if (!s) {
                    i = -i;
                    q = -q;
                }

                // Accumulation alternée pour casser la dépendance
                if (odd) {
                    accI[k][1] += i;
                    accQ[k][1] += q;
                } else {
                    accI[k][0] += i;
                    accQ[k][0] += q;
                }
            }
        }
    }

    // Finalisation: calcul de la puissance et émission
    OUTPUT_TILE: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS PIPELINE II=1
        int tau = tau_base + k;
        if (tau < NB_PHASES) {
            // Fusion des accumulateurs pair/impair
            acc_t accI_total = accI[k][0] + accI[k][1];
            acc_t accQ_total = accQ[k][0] + accQ[k][1];

            // Puissance = I² + Q²
            pow_pkt_t p;
            p.power  = (power_t)(accI_total * accI_total + accQ_total * accQ_total);
            p.tau    = tau;
            p.fd_idx = fd_idx;
            pow_out.write(p);
        }
    }
}

/**
 * Itère sur tous les tiles de phases pour un Doppler donné.
 */
static void process_tau_df_v9(
    const osc_t mixI_loc[N],
    const osc_t mixQ_loc[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<pow_pkt_t> &pow_out,
    int fd_idx
) {
#pragma HLS INLINE off

    // Parcours des phases par blocs de TAU_TILE
    TAU_TILE_LOOP: for (int tau_base = 0; tau_base < NB_PHASES; tau_base += TAU_TILE) {
    #pragma HLS LOOP_TRIPCOUNT min=128 max=128 avg=128
        process_tau_tile_v9(mixI_loc, mixQ_loc, prn_banks, tau_start_tbl, pow_out, fd_idx, tau_base);
    }
}

/**
 * Traitement complet pour une fréquence Doppler: mélange + corrélation.
 */
static void process_one_fd(
    const data_t signal_buf[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<pow_pkt_t> &pow_out,
    int fd_idx,
    int fd_step
) {
#pragma HLS INLINE off

    // Buffers locaux pour les données mélangées
    osc_t mixI_loc[N];
    osc_t mixQ_loc[N];
#pragma HLS ARRAY_PARTITION variable=mixI_loc cyclic factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=mixQ_loc cyclic factor=2 dim=1

    int fd_hz = FD_START + fd_idx * fd_step;

    doppler_mixer_to_mem(signal_buf, mixI_loc, mixQ_loc, fd_hz);
    process_tau_df_v9(mixI_loc, mixQ_loc, prn_banks, tau_start_tbl, pow_out, fd_idx);
}

// ============================================================
// ÉTAPE 5: GÉNÉRATION DE TOUTES LES PUISSANCES
// ============================================================

/**
 * Boucle sur toutes les fréquences Doppler et émet les puissances.
 */
static void produce_all_fd_powers(
    const data_t signal_buf[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    int nb_fd,
    int fd_step,
    hls::stream<pow_pkt_t> &pow_s
) {
#pragma HLS INLINE off

    FD_LOOP: for (int fd_idx = 0; fd_idx < nb_fd; fd_idx++) {
#pragma HLS LOOP_TRIPCOUNT min=1 max=21 avg=11
        process_one_fd(signal_buf, prn_banks, tau_start_tbl, pow_s, fd_idx, fd_step);
    }
}

// ============================================================
// ÉTAPE 6: RÉDUCTION (RECHERCHE DU MAXIMUM)
// ============================================================

/**
 * Consomme le stream de puissances, trouve le maximum global,
 * calcule la somme pour la moyenne, et écrit vers la sortie AXI-Stream.
 */
static void reduce_all_powers(
    hls::stream<pow_pkt_t> &pow_s,
    hls::stream<axis_t> &corr_out,
    int nb_fd,
    power_t &max_val,
    power_t &sum_corr,
    int &best_tau,
    int &best_fd_idx
) {
#pragma HLS INLINE off

    const int total = nb_fd * NB_PHASES;

    REDUCE_ALL: for (int i = 0; i < total; i++) {
#pragma HLS PIPELINE II=1
        pow_pkt_t p = pow_s.read();

        sum_corr += p.power;
        
        // Mise à jour du maximum
        if (p.power > max_val) {
            max_val     = p.power;
            best_tau    = p.tau;
            best_fd_idx = p.fd_idx;
        }

        // Émission vers PS via AXI-Stream
        axis_t out;
        out.data = (ap_int<32>)p.power;
        out.keep = -1;   // Tous les octets valides
        out.strb = -1;
        out.last = (i == total - 1) ? 1 : 0;  // TLAST sur le dernier
        corr_out.write(out);
    }
}

/**
 * Région dataflow principale: production et réduction en parallèle.
 */
static void run_fd_dataflow_region(
    const data_t signal_buf[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<axis_t> &corr_out,
    int nb_fd,
    int fd_step,
    power_t &max_val,
    power_t &sum_corr,
    int &best_tau,
    int &best_fd_idx
) {
#pragma HLS INLINE off

    // FIFO de liaison producteur-consommateur
    hls::stream<pow_pkt_t> pow_s("pow_s_top");
#pragma HLS STREAM variable=pow_s depth=64

#pragma HLS DATAFLOW
    produce_all_fd_powers(signal_buf, prn_banks, tau_start_tbl, nb_fd, fd_step, pow_s);
    reduce_all_powers(pow_s, corr_out, nb_fd, max_val, sum_corr, best_tau, best_fd_idx);
}

// ============================================================
// FONCTION TOP-LEVEL
// ============================================================

/**
 * Acquisition GPS série avec dataflow.
 * 
 * Algorithme:
 * 1. Capture signal GPS et code PRN depuis AXI-Stream
 * 2. Pour chaque fréquence Doppler [FD_START, FD_END]:
 *    a. Down-conversion via NCO (mélange I/Q)
 *    b. Corrélation avec PRN décalé pour toutes les phases tau
 * 3. Détection CFAR: pic > seuil = mean * (1 + K)
 * 4. Sortie: Doppler estimé, phase de code, flag détection
 */
void acquisition_serial(
    hls::stream<axis_t> &rx_real,
    hls::stream<axis_t> &prn_in,
    hls::stream<axis_t> &corr_out,
    int &doppler_out,
    int &codephase_out,
    int &sat_detected,
    int fd_step,
    int &max_power_out,
    int &mean_power_out,
    int &rx_count,
    int &prn_count,
    int &rx_last_seen,
    int &prn_last_seen,
    int &rx_last_pos,
    int &prn_last_pos
) {
    // Interfaces AXI
#pragma HLS INTERFACE axis port=rx_real
#pragma HLS INTERFACE axis port=prn_in
#pragma HLS INTERFACE axis port=corr_out
#pragma HLS INTERFACE s_axilite port=doppler_out   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=codephase_out bundle=CTRL
#pragma HLS INTERFACE s_axilite port=sat_detected  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=fd_step       bundle=CTRL
#pragma HLS INTERFACE s_axilite port=max_power_out  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=mean_power_out bundle=CTRL
#pragma HLS INTERFACE s_axilite port=return        bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_count bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_count bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_last_seen bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_last_seen bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_last_pos bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_last_pos bundle=CTRL

    // Buffers statiques (persistants entre appels, mappés en BRAM)
    static data_t signal_buf[N];                    // Signal GPS
    static ap_uint<1> prn_sign[2 * N];              // PRN étendu (×2 pour décalages)
    static ap_uint<1> prn_banks[TAU_TILE][2 * N];   // Copies PRN pour accès parallèle
    static int tau_start_tbl[NB_PHASES];            // Table des offsets de phase
#pragma HLS ARRAY_PARTITION variable=prn_sign cyclic factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=prn_banks complete dim=1

    // Validation du pas Doppler
    if (fd_step <= 0) {
        fd_step = 1;
    }

    int nb_fd = (FD_END - FD_START) / fd_step + 1;

    // Variables de recherche
    power_t max_val     = 0;
    power_t sum_corr    = 0;
    int best_tau        = 0;
    int best_fd_idx     = 0;

    // Reset debug TLAST
    rx_last_seen = 0;
    prn_last_seen = 0;
    rx_last_pos = -1;
    prn_last_pos = -1;

    // --- Phase 1: Capture des entrées ---
    load_inputs_once(
        rx_real, prn_in, signal_buf, prn_sign,
        rx_count, prn_count,
        rx_last_seen, prn_last_seen, rx_last_pos, prn_last_pos
    );
    
    // --- Phase 2: Construction des tables ---
    build_tau_start_table(tau_start_tbl);
    build_prn_banks(prn_sign, prn_banks);

    // --- Phase 3: Recherche grille Doppler×Phase ---
    run_fd_dataflow_region(
        signal_buf, prn_banks, tau_start_tbl, corr_out,
        nb_fd, fd_step,
        max_val, sum_corr, best_tau, best_fd_idx
    );

    // --- Phase 4: Détection CFAR ---
    int total_pts = nb_fd * NB_PHASES;
    power_t mean  = sum_corr / total_pts;

    // Seuil adaptatif: détection si max > mean × (1 + K)
    power_t threshold = mean * (power_t)(1.0 + SEUIL_K);
    sat_detected = (max_val > threshold) ? 1 : 0;

    // Sorties de monitoring
    max_power_out  = (int)max_val;
    mean_power_out = (int)mean;

    // Résultats d'acquisition
    doppler_out   = FD_START + best_fd_idx * fd_step;
    codephase_out = best_tau;
}
