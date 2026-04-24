#ifndef ACQ_SERIAL_DTFLOW_H
#define ACQ_SERIAL_DTFLOW_H

#include <hls_stream.h>
#include <ap_int.h>
#include <ap_fixed.h>
#include <hls_math.h>
#include <ap_axi_sdata.h>

// ============================================================
// PARAMÈTRES DE CONFIGURATION
// ============================================================

// Nombre d'échantillons par période de code GPS (1ms à ~12MHz)
#define N 11999 
#define FS 11999000.0
#define FREQUENCE_CENTRALE 3563000.0

// Seuils de détection CFAR
#define SEUIL_RATIO 2.5   // Ratio historique (non utilisé)
#define SEUIL_K    3.0    // CFAR: détection si max > mean * (1 + K)

// Plage de recherche Doppler (réduite en simulation pour accélérer)
#ifdef __SYNTHESIS__
    #define FD_START -10000
    #define FD_END   10000
    #define NB_PHASES 1023   // Nombre de phases de code à tester
#else
    #define FD_START -2500
    #define FD_END   -1000
    #define NB_PHASES 1023
#endif

// Configuration DDS (Direct Digital Synthesizer) pour NCO
#define DDS_LUT_BITS           10
#define DDS_LUT_SIZE           (1 << DDS_LUT_BITS)   // 1024 entrées LUT
#define DDS_PHASE_BITS         32                     // Résolution accumulateur de phase
#define FREQUENCE_CENTRALE_HZ  3563000
#define FS_INT                 11999000

#define NB_PHASES_TO_CHECK NB_PHASES
#define TWO_PI_DIV_FS  0.000000523641347317 

#ifndef ACQ_ENABLE_TLAST_DEBUG
#define ACQ_ENABLE_TLAST_DEBUG 0
#endif

// ============================================================
// DÉFINITION DES TYPES
// ============================================================

typedef ap_axis<32, 0, 0, 0> axis_t;  // AXI-Stream 32 bits

typedef ap_int<2> data_t;              // Signal GPS quantifié sur 2 bits

typedef ap_fixed<24, 8> osc_t;         // Sortie oscillateur NCO (16 bits frac)

typedef ap_fixed<32, 6> angle_t;       // Angle pour calculs trigonométriques
typedef ap_fixed<28, 2> trig_t;        // Résultat sin/cos

typedef ap_fixed<32, 24> acc_t;        // Accumulateur corrélation (8 bits frac)

typedef ap_fixed<64, 32> power_t;      // Puissance de corrélation (I² + Q²)

typedef ap_uint<32> phase_u_t;         // Accumulateur phase DDS (non signé)
typedef ap_int<32>  phase_t;           // Incrément phase DDS (signé)

// ============================================================
// STRUCTURES DE STREAMING DATAFLOW
// ============================================================

// Paquet signal reçu
struct rx_pkt_t {
    data_t     x;      // Échantillon
    ap_uint<1> last;   // Marqueur fin de trame
};

// Paquet après mélange NCO (down-conversion)
struct mix_pkt_t {
    osc_t      i;      // Voie en phase
    osc_t      q;      // Voie en quadrature
    ap_uint<1> last;
};

// Résultat de corrélation pour une cellule (tau, fd)
struct pow_pkt_t {
    power_t power;     // Puissance I² + Q²
    int     tau;       // Index phase de code [0, NB_PHASES-1]
    int     fd_idx;    // Index Doppler
};

// ============================================================
// PROTOTYPE FONCTION PRINCIPALE
// ============================================================

void acquisition_serial(
    hls::stream<axis_t> &rx_real,    // Entrée: signal GPS numérisé
    hls::stream<axis_t> &prn_in,     // Entrée: code PRN local
    hls::stream<axis_t> &corr_out,   // Sortie: grille de corrélation
    int &doppler_out,                // Sortie: Doppler du pic (Hz)
    int &codephase_out,              // Sortie: phase de code du pic
    int &sat_detected,               // Sortie: flag détection (0/1)
    int fd_step,                     // Entrée: pas Doppler (Hz)
    int &max_power_out,              // Sortie: puissance max (debug)
    int &mean_power_out,             // Sortie: puissance moyenne (debug)
    int &rx_count,                   // Debug: nb échantillons rx lus
    int &prn_count,                  // Debug: nb échantillons prn lus
    int &rx_last_seen,               // Debug: TLAST détecté sur rx
    int &prn_last_seen,              // Debug: TLAST détecté sur prn
    int &rx_last_pos,                // Debug: position TLAST rx
    int &prn_last_pos                // Debug: position TLAST prn
);

#endif
