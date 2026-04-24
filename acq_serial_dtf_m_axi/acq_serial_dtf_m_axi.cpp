/* #include "acq_serial_dtf_m_axi.h"
#include "dds_lut_rom.h"

static const int TAU_TILE = 8;

static inline ap_int<32> hz_to_phase_inc(int freq_hz) {
#pragma HLS INLINE
    const long long SCALE = (1LL << DDS_PHASE_BITS);
    long long num = (long long)freq_hz * SCALE;
    if (num >= 0) {
        num += (FS_INT / 2);
    } else {
        num -= (FS_INT / 2);
    }
    return (ap_int<32>)(num / (long long)FS_INT);
}

static void load_inputs_once(
    const mem_word_t *rx_real,
    const mem_word_t *prn_in,
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

    LOAD_INPUTS: for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        signal_buf[i] = (data_t)rx_real[i];
        prn_sign[i]   = (prn_in[i] > 0) ? (ap_uint<1>)1 : (ap_uint<1>)0;
    }

    rx_count      = N;
    prn_count     = N;
    rx_last_seen  = 0;
    prn_last_seen = 0;
    rx_last_pos   = -1;
    prn_last_pos  = -1;

    DUP_PRN: for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        prn_sign[i + N] = prn_sign[i];
    }
}

static void build_tau_start_table(int tau_start_tbl[NB_PHASES]) {
#pragma HLS INLINE off

    BUILD_TAU_START: for (int tau = 0; tau < NB_PHASES; tau++) {
#pragma HLS PIPELINE II=1
        tau_start_tbl[tau] = (int)(((long long)tau * (long long)N) / NB_PHASES);
    }
}

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

static void doppler_mixer_to_mem(
    const data_t signal_buf[N],
    osc_t mixI_loc[N],
    osc_t mixQ_loc[N],
    int fd_hz
) {
#pragma HLS INLINE off

    const ap_int<32> phase_inc = hz_to_phase_inc(-(FREQUENCE_CENTRALE_HZ + fd_hz));
    ap_uint<32> phase_acc = 0;

    MIX_ALL: for (int n = 0; n < N; n++) {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=mixI_loc inter false
#pragma HLS DEPENDENCE variable=mixQ_loc inter false
        ap_uint<DDS_LUT_BITS> lut_idx =
            phase_acc.range(DDS_PHASE_BITS - 1, DDS_PHASE_BITS - DDS_LUT_BITS);

        osc_t c = DDS_COS_LUT[(int)lut_idx];
        osc_t s = DDS_SIN_LUT[(int)lut_idx];
        osc_t x = (osc_t)signal_buf[n];

        mixI_loc[n] = x * c;
        mixQ_loc[n] = -(x * s);
        phase_acc   = (ap_uint<32>)(phase_acc + (ap_uint<32>)phase_inc);
    }
}

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

    acc_t accI[TAU_TILE][2];
    acc_t accQ[TAU_TILE][2];
#pragma HLS ARRAY_PARTITION variable=accI complete dim=1
#pragma HLS ARRAY_PARTITION variable=accI complete dim=2
#pragma HLS ARRAY_PARTITION variable=accQ complete dim=1
#pragma HLS ARRAY_PARTITION variable=accQ complete dim=2

    INIT_TILE: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS UNROLL
        accI[k][0] = 0;
        accI[k][1] = 0;
        accQ[k][0] = 0;
        accQ[k][1] = 0;
    }

    SAMPLE_LOOP: for (int n = 0; n < N; n++) {
#pragma HLS PIPELINE II=1
        osc_t base_i = mixI_loc[n];
        osc_t base_q = mixQ_loc[n];
        bool odd = (n & 1) != 0;

        TAU_TILE_LOOP: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS UNROLL
            int tau = tau_base + k;
            if (tau < NB_PHASES) {
                osc_t i = base_i;
                osc_t q = base_q;
                ap_uint<1> s = prn_banks[k][tau_start_tbl[tau] + n];

                if (!s) {
                    i = -i;
                    q = -q;
                }

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

    OUTPUT_TILE: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS PIPELINE II=1
        int tau = tau_base + k;
        if (tau < NB_PHASES) {
            acc_t accI_total = accI[k][0] + accI[k][1];
            acc_t accQ_total = accQ[k][0] + accQ[k][1];

            pow_pkt_t p;
            p.power  = (power_t)(accI_total * accI_total + accQ_total * accQ_total);
            p.tau    = tau;
            p.fd_idx = fd_idx;
            pow_out.write(p);
        }
    }
}

static void process_tau_df_v9(
    const osc_t mixI_loc[N],
    const osc_t mixQ_loc[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<pow_pkt_t> &pow_out,
    int fd_idx
) {
#pragma HLS INLINE off

    TAU_TILE_LOOP: for (int tau_base = 0; tau_base < NB_PHASES; tau_base += TAU_TILE) {
#pragma HLS LOOP_TRIPCOUNT min=128 max=128 avg=128
        process_tau_tile_v9(mixI_loc, mixQ_loc, prn_banks, tau_start_tbl, pow_out, fd_idx, tau_base);
    }
}

static void process_one_fd(
    const data_t signal_buf[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<pow_pkt_t> &pow_out,
    int fd_idx,
    int fd_step
) {
#pragma HLS INLINE off

    osc_t mixI_loc[N];
    osc_t mixQ_loc[N];
#pragma HLS ARRAY_PARTITION variable=mixI_loc cyclic factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=mixQ_loc cyclic factor=2 dim=1

    int fd_hz = FD_START + fd_idx * fd_step;

    doppler_mixer_to_mem(signal_buf, mixI_loc, mixQ_loc, fd_hz);
    process_tau_df_v9(mixI_loc, mixQ_loc, prn_banks, tau_start_tbl, pow_out, fd_idx);
}

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

static void reduce_all_powers(
    hls::stream<pow_pkt_t> &pow_s,
    mem_word_t *corr_out,
    int nb_fd,
    power_t &max_val,
    power_t &sum_corr,
    int &best_tau,
    int &best_fd_idx,
    int &corr_count
) {
#pragma HLS INLINE off

    const int total = nb_fd * NB_PHASES;

    REDUCE_ALL: for (int i = 0; i < total; i++) {
#pragma HLS PIPELINE II=1
        pow_pkt_t p = pow_s.read();

        sum_corr += p.power;
        if (p.power > max_val) {
            max_val     = p.power;
            best_tau    = p.tau;
            best_fd_idx = p.fd_idx;
        }

        corr_out[i] = (mem_word_t)p.power;
    }

    corr_count = total;
}

static void run_fd_dataflow_region(
    const data_t signal_buf[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    mem_word_t *corr_out,
    int nb_fd,
    int fd_step,
    power_t &max_val,
    power_t &sum_corr,
    int &best_tau,
    int &best_fd_idx,
    int &corr_count
) {
#pragma HLS INLINE off

    hls::stream<pow_pkt_t> pow_s("pow_s_top");
#pragma HLS STREAM variable=pow_s depth=64

#pragma HLS DATAFLOW
    produce_all_fd_powers(signal_buf, prn_banks, tau_start_tbl, nb_fd, fd_step, pow_s);
    reduce_all_powers(pow_s, corr_out, nb_fd, max_val, sum_corr, best_tau, best_fd_idx, corr_count);
}

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
) {
#pragma HLS INTERFACE m_axi port=rx_real  offset=slave bundle=GMEM0 depth=N
#pragma HLS INTERFACE m_axi port=prn_in   offset=slave bundle=GMEM1 depth=N
#pragma HLS INTERFACE m_axi port=corr_out offset=slave bundle=GMEM2

#pragma HLS INTERFACE s_axilite port=rx_real        bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_in         bundle=CTRL
#pragma HLS INTERFACE s_axilite port=corr_out       bundle=CTRL
#pragma HLS INTERFACE s_axilite port=corr_count     bundle=CTRL
#pragma HLS INTERFACE s_axilite port=doppler_out    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=codephase_out  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=sat_detected   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=fd_step        bundle=CTRL
#pragma HLS INTERFACE s_axilite port=max_power_out  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=mean_power_out bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_count       bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_count      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_last_seen   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_last_seen  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_last_pos    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_last_pos   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=return         bundle=CTRL

    static data_t signal_buf[N];
    static ap_uint<1> prn_sign[2 * N];
    static ap_uint<1> prn_banks[TAU_TILE][2 * N];
    static int tau_start_tbl[NB_PHASES];
#pragma HLS ARRAY_PARTITION variable=prn_sign cyclic factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=prn_banks complete dim=1

    if (fd_step <= 0) {
        fd_step = 1;
    }

    int nb_fd = (FD_END - FD_START) / fd_step + 1;

    power_t max_val  = 0;
    power_t sum_corr = 0;
    int best_tau     = 0;
    int best_fd_idx  = 0;

    corr_count = 0;

    load_inputs_once(
        rx_real,
        prn_in,
        signal_buf,
        prn_sign,
        rx_count,
        prn_count,
        rx_last_seen,
        prn_last_seen,
        rx_last_pos,
        prn_last_pos
    );

    build_tau_start_table(tau_start_tbl);
    build_prn_banks(prn_sign, prn_banks);

    run_fd_dataflow_region(
        signal_buf,
        prn_banks,
        tau_start_tbl,
        corr_out,
        nb_fd,
        fd_step,
        max_val,
        sum_corr,
        best_tau,
        best_fd_idx,
        corr_count
    );

    int total_pts = nb_fd * NB_PHASES;
    power_t mean  = sum_corr / total_pts;

    power_t threshold = mean * (power_t)(1.0 + SEUIL_K);
    sat_detected = (max_val > threshold) ? 1 : 0;

    max_power_out  = (int)max_val;
    mean_power_out = (int)mean;
    doppler_out    = FD_START + best_fd_idx * fd_step;
    codephase_out  = best_tau;
} */

















#include "acq_serial_dtf_m_axi.h"
#include "dds_lut_rom.h"

static const int TAU_TILE = 8;

// ========================================================
// INPUT DATAFLOW (AXI FULL) : lecture memoire -> streams -> buffers
// ========================================================
static void read_inputs_to_stream_m_axi(
    const mem_word_t *rx_real,
    const mem_word_t *prn_in,
    hls::stream<data_t> &sig_s,
    hls::stream<ap_uint<1> > &prn_s
) {
#pragma HLS INLINE off

    READ_INPUTS: for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        mem_word_t rxw  = rx_real[i];
        mem_word_t prnw = prn_in[i];

        sig_s.write((data_t)rxw);
        prn_s.write((prnw > 0) ? (ap_uint<1>)1 : (ap_uint<1>)0);
    }
}

static void store_stream_to_mem(
    hls::stream<data_t> &sig_s,
    hls::stream<ap_uint<1> > &prn_s,
    data_t signal_buf[N],
    ap_uint<1> prn_sign[2 * N]
) {
#pragma HLS INLINE off

    STORE_INPUTS: for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        signal_buf[i] = sig_s.read();
        prn_sign[i]   = prn_s.read();
    }
}

static void capture_inputs_to_mem_m_axi(
    const mem_word_t *rx_real,
    const mem_word_t *prn_in,
    data_t signal_buf[N],
    ap_uint<1> prn_sign[2 * N]
) {
#pragma HLS INLINE off

    hls::stream<data_t> sig_s("sig_s");
    hls::stream<ap_uint<1> > prn_s("prn_s");
#pragma HLS STREAM variable=sig_s depth=64
#pragma HLS STREAM variable=prn_s depth=64

#pragma HLS DATAFLOW
    read_inputs_to_stream_m_axi(rx_real, prn_in, sig_s, prn_s);
    store_stream_to_mem(sig_s, prn_s, signal_buf, prn_sign);
}

static inline ap_int<32> hz_to_phase_inc(int freq_hz) {
#pragma HLS INLINE
    const long long SCALE = (1LL << DDS_PHASE_BITS);
    long long num = (long long)freq_hz * SCALE;
    if (num >= 0) {
        num += (FS_INT / 2);
    } else {
        num -= (FS_INT / 2);
    }
    return (ap_int<32>)(num / (long long)FS_INT);
}

static void load_inputs_once(
    const mem_word_t *rx_real,
    const mem_word_t *prn_in,
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

    rx_count = 0;
    prn_count = 0;

    capture_inputs_to_mem_m_axi(
        rx_real,
        prn_in,
        signal_buf,
        prn_sign
    );

    rx_count      = N;
    prn_count     = N;
    rx_last_seen  = 0;   // pas de TLAST en AXI Full
    prn_last_seen = 0;   // pas de TLAST en AXI Full
    rx_last_pos   = -1;
    prn_last_pos  = -1;

    DUP_PRN: for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1
        prn_sign[i + N] = prn_sign[i];
    }
}

static void build_tau_start_table(int tau_start_tbl[NB_PHASES]) {
#pragma HLS INLINE off

    BUILD_TAU_START: for (int tau = 0; tau < NB_PHASES; tau++) {
#pragma HLS PIPELINE II=1
        tau_start_tbl[tau] = (int)(((long long)tau * (long long)N) / NB_PHASES);
    }
}

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

static void doppler_mixer_to_mem(
    const data_t signal_buf[N],
    osc_t mixI_loc[N],
    osc_t mixQ_loc[N],
    int fd_hz
) {
#pragma HLS INLINE off

    const ap_int<32> phase_inc = hz_to_phase_inc(-(FREQUENCE_CENTRALE_HZ + fd_hz));
    ap_uint<32> phase_acc = 0;

    MIX_ALL: for (int n = 0; n < N; n++) {
#pragma HLS PIPELINE II=1
#pragma HLS DEPENDENCE variable=mixI_loc inter false
#pragma HLS DEPENDENCE variable=mixQ_loc inter false
        ap_uint<DDS_LUT_BITS> lut_idx =
            phase_acc.range(DDS_PHASE_BITS - 1, DDS_PHASE_BITS - DDS_LUT_BITS);

        osc_t c = DDS_COS_LUT[(int)lut_idx];
        osc_t s = DDS_SIN_LUT[(int)lut_idx];
        osc_t x = (osc_t)signal_buf[n];

        mixI_loc[n] = x * c;
        mixQ_loc[n] = -(x * s);
        phase_acc   = (ap_uint<32>)(phase_acc + (ap_uint<32>)phase_inc);
    }
}

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

    acc_t accI[TAU_TILE][2];
    acc_t accQ[TAU_TILE][2];
#pragma HLS ARRAY_PARTITION variable=accI complete dim=1
#pragma HLS ARRAY_PARTITION variable=accI complete dim=2
#pragma HLS ARRAY_PARTITION variable=accQ complete dim=1
#pragma HLS ARRAY_PARTITION variable=accQ complete dim=2

    INIT_TILE: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS UNROLL
        accI[k][0] = 0;
        accI[k][1] = 0;
        accQ[k][0] = 0;
        accQ[k][1] = 0;
    }

    SAMPLE_LOOP: for (int n = 0; n < N; n++) {
#pragma HLS PIPELINE II=1
        osc_t base_i = mixI_loc[n];
        osc_t base_q = mixQ_loc[n];
        bool odd = (n & 1) != 0;

        TAU_TILE_LOOP: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS UNROLL
            int tau = tau_base + k;
            if (tau < NB_PHASES) {
                osc_t i = base_i;
                osc_t q = base_q;
                ap_uint<1> s = prn_banks[k][tau_start_tbl[tau] + n];

                if (!s) {
                    i = -i;
                    q = -q;
                }

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

    OUTPUT_TILE: for (int k = 0; k < TAU_TILE; k++) {
#pragma HLS PIPELINE II=1
        int tau = tau_base + k;
        if (tau < NB_PHASES) {
            acc_t accI_total = accI[k][0] + accI[k][1];
            acc_t accQ_total = accQ[k][0] + accQ[k][1];

            pow_pkt_t p;
            p.power  = (power_t)(accI_total * accI_total + accQ_total * accQ_total);
            p.tau    = tau;
            p.fd_idx = fd_idx;
            pow_out.write(p);
        }
    }
}

static void process_tau_df_v9(
    const osc_t mixI_loc[N],
    const osc_t mixQ_loc[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<pow_pkt_t> &pow_out,
    int fd_idx
) {
#pragma HLS INLINE off

    TAU_TILE_LOOP: for (int tau_base = 0; tau_base < NB_PHASES; tau_base += TAU_TILE) {
#pragma HLS LOOP_TRIPCOUNT min=128 max=128 avg=128
        process_tau_tile_v9(mixI_loc, mixQ_loc, prn_banks, tau_start_tbl, pow_out, fd_idx, tau_base);
    }
}

static void process_one_fd(
    const data_t signal_buf[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    hls::stream<pow_pkt_t> &pow_out,
    int fd_idx,
    int fd_step
) {
#pragma HLS INLINE off

    osc_t mixI_loc[N];
    osc_t mixQ_loc[N];
#pragma HLS ARRAY_PARTITION variable=mixI_loc cyclic factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=mixQ_loc cyclic factor=2 dim=1

    int fd_hz = FD_START + fd_idx * fd_step;

    doppler_mixer_to_mem(signal_buf, mixI_loc, mixQ_loc, fd_hz);
    process_tau_df_v9(mixI_loc, mixQ_loc, prn_banks, tau_start_tbl, pow_out, fd_idx);
}

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

static void reduce_all_powers(
    hls::stream<pow_pkt_t> &pow_s,
    mem_word_t *corr_out,
    int nb_fd,
    power_t &max_val,
    power_t &sum_corr,
    int &best_tau,
    int &best_fd_idx,
    int &corr_count
) {
#pragma HLS INLINE off

    const int total = nb_fd * NB_PHASES;

    REDUCE_ALL: for (int i = 0; i < total; i++) {
#pragma HLS PIPELINE II=1
        pow_pkt_t p = pow_s.read();

        sum_corr += p.power;
        if (p.power > max_val) {
            max_val     = p.power;
            best_tau    = p.tau;
            best_fd_idx = p.fd_idx;
        }

        corr_out[i] = (mem_word_t)p.power;
    }

    corr_count = total;
}

static void run_fd_dataflow_region(
    const data_t signal_buf[N],
    const ap_uint<1> prn_banks[TAU_TILE][2 * N],
    const int tau_start_tbl[NB_PHASES],
    mem_word_t *corr_out,
    int nb_fd,
    int fd_step,
    power_t &max_val,
    power_t &sum_corr,
    int &best_tau,
    int &best_fd_idx,
    int &corr_count
) {
#pragma HLS INLINE off

    hls::stream<pow_pkt_t> pow_s("pow_s_top");
#pragma HLS STREAM variable=pow_s depth=64

#pragma HLS DATAFLOW
    produce_all_fd_powers(signal_buf, prn_banks, tau_start_tbl, nb_fd, fd_step, pow_s);
    reduce_all_powers(pow_s, corr_out, nb_fd, max_val, sum_corr, best_tau, best_fd_idx, corr_count);
}

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
) {
#pragma HLS INTERFACE m_axi port=rx_real  offset=slave bundle=GMEM0 depth=N
#pragma HLS INTERFACE m_axi port=prn_in   offset=slave bundle=GMEM1 depth=N
#pragma HLS INTERFACE m_axi port=corr_out offset=slave bundle=GMEM2

#pragma HLS INTERFACE s_axilite port=rx_real        bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_in         bundle=CTRL
#pragma HLS INTERFACE s_axilite port=corr_out       bundle=CTRL
#pragma HLS INTERFACE s_axilite port=corr_count     bundle=CTRL
#pragma HLS INTERFACE s_axilite port=doppler_out    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=codephase_out  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=sat_detected   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=fd_step        bundle=CTRL
#pragma HLS INTERFACE s_axilite port=max_power_out  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=mean_power_out bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_count       bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_count      bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_last_seen   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_last_seen  bundle=CTRL
#pragma HLS INTERFACE s_axilite port=rx_last_pos    bundle=CTRL
#pragma HLS INTERFACE s_axilite port=prn_last_pos   bundle=CTRL
#pragma HLS INTERFACE s_axilite port=return         bundle=CTRL

    static data_t signal_buf[N];
    static ap_uint<1> prn_sign[2 * N];
    static ap_uint<1> prn_banks[TAU_TILE][2 * N];
    static int tau_start_tbl[NB_PHASES];
#pragma HLS ARRAY_PARTITION variable=prn_sign cyclic factor=2 dim=1
#pragma HLS ARRAY_PARTITION variable=prn_banks complete dim=1

    if (fd_step <= 0) {
        fd_step = 1;
    }

    int nb_fd = (FD_END - FD_START) / fd_step + 1;

    power_t max_val  = 0;
    power_t sum_corr = 0;
    int best_tau     = 0;
    int best_fd_idx  = 0;

    corr_count = 0;

    load_inputs_once(
        rx_real,
        prn_in,
        signal_buf,
        prn_sign,
        rx_count,
        prn_count,
        rx_last_seen,
        prn_last_seen,
        rx_last_pos,
        prn_last_pos
    );

    build_tau_start_table(tau_start_tbl);
    build_prn_banks(prn_sign, prn_banks);

    run_fd_dataflow_region(
        signal_buf,
        prn_banks,
        tau_start_tbl,
        corr_out,
        nb_fd,
        fd_step,
        max_val,
        sum_corr,
        best_tau,
        best_fd_idx,
        corr_count
    );

    int total_pts = nb_fd * NB_PHASES;
    power_t mean  = sum_corr / total_pts;

    power_t threshold = mean * (power_t)(1.0 + SEUIL_K);
    sat_detected = (max_val > threshold) ? 1 : 0;

    max_power_out  = (int)max_val;
    mean_power_out = (int)mean;
    doppler_out    = FD_START + best_fd_idx * fd_step;
    codephase_out  = best_tau;
}