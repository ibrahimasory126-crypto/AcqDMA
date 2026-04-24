#ifndef ACQ_SERIAL_DTFLOW_H
#define ACQ_SERIAL_DTFLOW_H

#include <hls_stream.h>
#include <ap_int.h>
#include <ap_fixed.h>
#include <hls_math.h>

// --- Paramètres FIXES ---
#define N 11999
#define FS 11999000.0
#define FREQUENCE_CENTRALE 3563000.0
#define SEUIL_RATIO 2.5
#define SEUIL_K 3.0

// --- Paramètres ADAPTATIFS (CSIM vs SYNTHESE) ---
#ifdef __SYNTHESIS__
    #define FD_START -10000
    #define FD_END   10000
    #define NB_PHASES 1023
#else
    #define FD_START -2500
    #define FD_END   -1000
    #define NB_PHASES 1023
#endif

// --- Paramètres DDS LUT ---
#define DDS_LUT_BITS           10
#define DDS_LUT_SIZE           (1 << DDS_LUT_BITS)
#define DDS_PHASE_BITS         32
#define FREQUENCE_CENTRALE_HZ  3563000
#define FS_INT                 11999000

#define NB_PHASES_TO_CHECK NB_PHASES
#define TWO_PI_DIV_FS 0.000000523641347317

#ifndef ACQ_ENABLE_TLAST_DEBUG
#define ACQ_ENABLE_TLAST_DEBUG 0
#endif

typedef ap_int<32> mem_word_t;

// Types de données
typedef ap_int<2> data_t;
typedef ap_fixed<24, 8> osc_t;
typedef ap_fixed<32, 6> angle_t;
typedef ap_fixed<28, 2> trig_t;
typedef ap_fixed<32, 24> acc_t;
typedef ap_fixed<64, 32> power_t;
typedef ap_uint<32> phase_u_t;
typedef ap_int<32>  phase_t;

struct rx_pkt_t {
    data_t     x;
    ap_uint<1> last;
};

struct mix_pkt_t {
    osc_t      i;
    osc_t      q;
    ap_uint<1> last;
};

struct pow_pkt_t {
    power_t power;
    int     tau;
    int     fd_idx;
};

void acquisition_serial_m_axi(
    const mem_word_t *rx_real,
    const mem_word_t *prn_in,
    mem_word_t *corr_out,
    int &corr_count,
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
);

#endif