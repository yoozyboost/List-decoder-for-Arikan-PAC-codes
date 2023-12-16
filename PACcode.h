#ifndef PAC_PACCODE_H
#define PAC_PACCODE_H


#include <cstdint>
#include <vector>
#include <math.h>
#include <stack>          // std::stack
#include <cstdio>
#include <algorithm>

class PACcode{

    public:

        PACcode(uint16_t block_length,uint16_t info_length,uint16_t construction_mode,double design_snr,int crc_mode,uint16_t S,std::vector<uint16_t> segmented_bits_length,std::vector<uint16_t> segmented_crc_length,std::vector<uint16_t> segmented_list_size):
            _block_length(block_length), _info_length(info_length),_design_snr(design_snr),_crc_mode(crc_mode),_construction_mode(construction_mode), _S(S), _segmented_bits_length(segmented_bits_length),_segmented_crc_length(segmented_crc_length),_segmented_list_size(segmented_list_size)
        {
            _n = (uint16_t)log2(_block_length);
            _bit_rev_order.resize(_block_length);
            _frozen_bits.resize(block_length);
            _channel_order_descending.resize(_block_length);

            //generator polynomial for convolutional encoder
            _m = 7;
            convolutional_polynomial = {1,1,0,1,1,0,1};
        
            create_bit_rev_order();
            initialize_frozen_bits();

            _max_list_size = *max_element(_segmented_list_size.begin(),_segmented_list_size.end());

            //initialize _shift_register of 2D array of size L*(_m)
            _shift_register.resize(_max_list_size);
            for (int i = 0; i < _max_list_size; i++)
            {
                _shift_register.at(i).resize(_m);
            }

        }

        //general functions for simulation
        std::vector<uint8_t> dataGenerator(uint16_t data_length);
        std::vector<int8_t> BPSKmodulation(std::vector<uint8_t> data_to_modulate);
        std::vector<double> AWGNchannel(std::vector<int8_t> send_data,double dispersion);
        std::vector<uint8_t> BPSKdemodulation(std::vector<double> data_to_demodulate);
        uint16_t errorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2);
        uint16_t blockErrorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2);


        //PAC codes functions
        std::vector<uint8_t> rate_profile(std::vector<uint8_t> info_bits);
        std::vector<uint8_t> convolutional_encoder(std::vector<uint8_t> profiled_bits);

        //polar codes functions
        std::vector<uint8_t> polar_encoder(std::vector<uint8_t> info_bits);
        std::vector<uint8_t> decode_scascl_llr(std::vector<double> llr, uint16_t list_size);

    private:
        uint16_t _m;
        uint16_t _n;
        uint16_t _S;
        uint16_t _block_length;
        uint16_t _info_length;
        uint16_t _crc_size;
        double _design_snr;

        uint16_t _current_segment;
        uint16_t _segmented_block_first_index;
        uint16_t _segmented_block_last_index;
        
        uint16_t _max_list_size;
        uint16_t _previous_list_size;
        uint16_t _current_list_size; 

        uint16_t _construction_mode;
        uint16_t _crc_mode;

        std::vector<uint8_t> _crc_array;

        std::vector<uint16_t> _segmented_bits_length;
        std::vector<uint16_t> _segmented_crc_length;
        std::vector<uint16_t> _segmented_list_size;

        std::vector<uint16_t> _bit_rev_order;
        std::vector<uint16_t> _channel_order_descending;
        std::vector<uint8_t> _frozen_bits;

        std::vector<uint8_t> _partial_unfrozen_bits;

        std::vector<uint8_t> _sending_bits;
        std::vector<uint8_t> _received_bits;

        std::vector<uint8_t> _priority_info_bits;
        std::vector<uint8_t> _priority_decoded_info_bits;

        std::vector<std::vector<uint8_t>> _segmented_info_bits;
        std::vector<std::vector<uint8_t>> _segmented_decoded_info_bits;

        std::vector<uint8_t> _regular_info_bits;
        std::vector<uint8_t> _regular_decoded_info_bits;

        std::vector<std::vector<uint8_t>> _shift_register;
        std::vector<uint8_t> convolutional_polynomial;     

        void crc_type(uint16_t crc_size);
        void initialize_frozen_bits();
        void create_bit_rev_order();

        std::vector<uint8_t> decode_scl();
        std::vector<uint8_t> decode_scascl();
        bool _llr_based_computation;

        std::vector<std::vector<double *>> _arrayPointer_LLR;
        std::vector<double> _pathMetric_LLR;

        uint16_t _list_size;

        std::stack<uint16_t> _inactivePathIndices;
        std::vector<uint16_t > _activePath;
        std::vector<std::vector<double *>> _arrayPointer_P;
        std::vector<std::vector<uint8_t *>> _arrayPointer_C;
        std::vector<uint8_t *> _arrayPointer_Info;
        std::vector<std::vector<uint16_t>> _pathIndexToArrayIndex;
        std::vector<std::stack<uint16_t>> _inactiveArrayIndices;
        std::vector<std::vector<uint16_t>> _arrayReferenceCount;

        std::vector<std::vector<uint8_t>> _priorityPaths;

        std::vector<uint8_t> _crc_bits;
        std::vector<int> _list_correct_number;


        void initializeDataStructures();
        uint16_t assignInitialPath();
        uint16_t clonePath(uint16_t);
        void killPath(uint16_t l);

        double * getArrayPointer_P(uint16_t lambda, uint16_t  l);
        double * getArrayPointer_LLR(uint16_t lambda, uint16_t  l);
        uint8_t * getArrayPointer_C(uint16_t lambda, uint16_t  l);

        void recursivelyCalcP(uint16_t lambda, uint16_t phi);
        void recursivelyCalcLLR(uint16_t lambda, uint16_t phi);
        void recursivelyUpdateC(uint16_t lambda, uint16_t phi);

        void continuePaths_FrozenBit(uint16_t phi);
        void continuePaths_UnfrozenBit(uint16_t phi);

        uint16_t findMostProbablePath(bool check_crc);
        
        bool crc_check(uint8_t * info_bits_padded); 
        uint8_t segment_crc_check(std::vector<uint8_t> info_bit_padded);

        void change_segment();

        double in_out(double x);
        double calcGamma(double u);
        double calcXi(double gamma);
        double calcXiInverse(double z);
        double bisection(double y,double a,double b,double eps);
        double ff(double x,double y);
        double CALC_CAP(double xi);
        
};

#endif //PAC_PACCODE_H