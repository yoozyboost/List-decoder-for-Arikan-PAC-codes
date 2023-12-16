#include "PACcode.h"
#include <iostream>
#include <cmath>       /* log */
#include <sstream>      // std::stringstream
#include <fstream>
#include <iomanip>      // std::setprecision
#include <random>
#include <algorithm>
#include <chrono>
#include <bitset>
#include <map>

std::vector<uint8_t> PACcode::dataGenerator(uint16_t data_length){

    std::vector<uint8_t> generated_data(data_length,0);

    for(uint16_t i=0;i<data_length;i++){
        generated_data.at(i) = (rand()%2);
    }

    return generated_data;
}


std::vector<double> PACcode::AWGNchannel(std::vector<int8_t> send_data, double dispersion){

    uint16_t data_size = send_data.size();
    std::vector<double> sent_data(data_size);

    double u1,u2;
    double v1,v2;

    double sqrt2inv = 1.0/sqrt(2.0);

    for(uint16_t i = 0; i < data_size; i++ ){
		// generating random values
		do{
			u1 = (double) rand()/RAND_MAX;
		} while( u1 == 0 );
		u2 = (double) rand()/RAND_MAX;

		v1 = sqrt( -2.0 * log( u1 ) );
		v2 = 2.0 * M_PI * u2;

		sent_data.at(i) = send_data.at(i) +  v1 * dispersion * cos(v2) * sqrt2inv;
    }

    return sent_data;

}

std::vector<int8_t> PACcode::BPSKmodulation(std::vector<uint8_t> data_to_modulate){
    uint16_t data_size = data_to_modulate.size();
    std::vector<int8_t> modulated_data(data_size);

    for(int16_t i=0;i<data_size;i++){
        modulated_data.at(i) = 1 - 2 * data_to_modulate.at(i);
        // modulated_data.at(i+1) = 1 - 2 * data_to_modulate.at(i+1);
    }

    return modulated_data;

}

std::vector<uint8_t> PACcode::BPSKdemodulation(std::vector<double> data_to_demodulate){
    uint16_t data_size = data_to_demodulate.size();
    std::vector<uint8_t> received_data(data_size);

    for(int16_t i=0;i<data_size;i++){
        if(data_to_demodulate.at(i) > 0.0){
            received_data.at(i) = 0;
        }
        else{
            received_data.at(i) = 1;
        }
    }

    return received_data;

}

uint16_t PACcode::errorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2){
    uint16_t data_size = data1.size();
    uint16_t error_count = 0;
    for(uint16_t i=0;i<data_size;i++){
        if(data1.at(i) != data2.at(i)){
            error_count++;
        }
    }
    return error_count;
}

uint16_t PACcode::blockErrorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2){
    uint16_t data_size = data1.size();
    uint16_t error_count = 0;
    for(uint16_t i=0;i<data_size;i++){
        if(data1.at(i) != data2.at(i)){
            error_count = 1;
            break;
        }
    }
    return error_count;
}

void PACcode::crc_type(uint16_t crc_size){

    std::map <uint16_t, std::vector<uint8_t>> crc_map = {
    {64, {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1}},
    {48, {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,1,1}},
    {32, {1,0,0,1,0,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,1,1,1}},
    {31, {1,0,0,0,1,1,0,1,1,1,0,0,1,0,1,0,1,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1}},
    {30, {1,0,0,0,1,0,1,0,0,1,1,1,0,1,1,1,1,1,0,0,0,1,1,0,1,0,1,1,0,0,1}},
    {29, {1,0,1,0,1,1,1,1,0,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,1}},

    //0x91dc1e3
    {28, {1,0,0,1,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1}},
    {27, {1,1,0,0,1,0,1,1,0,1,1,1,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,1}},

    //0x33c19ef
    {26, {1,1,0,0,1,1,1,1,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,1,1}},
    {25, {1,0,0,1,0,1,0,0,1,0,0,0,1,1,1,0,0,1,1,0,0,1,1,1,1,1}},


    //0x8f90e3
    // {24, {1,0,0,0,1,1,1,1,1,0,0,1,0,0,0,0,1,1,1,0,0,0,1,1,1}},
    //0xb73e91
    // {24, {1,0,1,1,0,1,1,1,0,0,1,1,1,1,1,0,1,0,0,1,0,0,0,1,1}},

    //0xC3267D
    {24, {1,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1}},

    //0x540df0
    {23, {1,0,1,0,1,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,0,0,0,1}},
    //0x4b79d1
    // {23, {1,0,0,1,0,1,1,0,1,1,1,1,0,0,1,1,1,0,1,0,0,0,1,1}},

    //0x25d467
    // {22, {1,0,0,1,0,1,1,1,0,1,0,1,0,0,0,1,1,0,0,1,1,1,1}},

    //0x308fd3
    {22, {1,1,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,1,0,0,1,1,1}},

    //0x165751
    // {21, {1,0,1,1,0,0,1,0,1,0,1,1,1,0,1,0,0,0,1,0,1,1}},

    //0x1707ea
    {21, {1,0,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,1,0,1,0,1}},

    //0x9d587
    // {20, {1,0,0,1,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,1,1}},

    //0xb5827
    {20, {1,0,1,1,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,1}},

    {19, {1,0,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,0,1,1}},

    //0x2235b
    {18, {1,0,0,0,1,0,0,0,1,1,0,1,0,1,1,0,1,1,1}},
    {17, {1,0,1,1,1,0,1,1,0,1,0,1,0,0,1,1,1,1}},

    // 
    {16, {1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1}},

    //0xd157
    // {16, {1,1,0,1,0,0,0,1,0,1,1,1,0,1,0,1,1}},

    {15, {1,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1}},

    //0x372b
    {14, {1,1,0,1,1,1,0,0,1,0,1,0,1,1,1}},

    {13, {1,0,0,0,0,1,0,1,1,0,1,1,1,1}},


    // {12, {1,0,1,0,0,1,0,0,1,1,1,1,1}},

    //0xc07
    // {12, {1,1,0,0,0,0,0,0,0,1,1,1,1}},
    //0xcab
    {12, {1,1,0,0,1,0,1,0,1,0,1,1,1}},


    // {11, {1,0,0,1,1,1,1,0,1,0,1,1}},

    //0x49f
    {11, {1,0,0,1,0,0,1,1,1,1,1,1}},

    //0x327
    {10, {1,1,0,0,1,0,0,1,1,1,1}},

    //0x319
    // {10, {1,1,0,0,0,1,1,0,0,1,1}},

    //0x2b9
    // {10, {1,0,1,0,1,1,1,0,0,1,1}},

    //0x247
    // {10, {1,0,0,1,0,0,0,1,1,1,1}},

    //0x28e
    // {10, {1,0,1,0,0,0,1,1,1,0,1}},

    //0x29b
    //  {10, {1,0,1,0,0,1,1,0,1,1,1}},
    

    //0x233
    // {10, {1,0,0,0,1,1,0,0,1,1,1}},

    //0x3D9
    // {10, {1,1,1,1,0,1,1,0,0,1,1}},

    //0x34C
    // {10, {1,1,0,1,0,0,1,1,0,0,1}},

    //0x13c
    {9,  {1,0,0,1,1,1,1,0,0,1}},

    //0x185
    // {9,  {1,1,0,0,0,0,1,0,1,1}},

    // 0x119
    // {9,  {1,0,0,0,1,1,0,0,1,1}},

    // 0x14b
    // {9,  {1,0,1,0,0,1,0,1,1,1}},

    // 0x167
    // {9,  {1,0,1,1,0,0,1,1,1,1}},

    
    //0xA6
    {8,  {1,0,1,0,0,1,1,0,1}},

    //0x97
    // {8,  {1,0,0,1,0,1,1,1,1}},
    
    {7,  {1,1,1,0,0,1,0,1}},

    //0x21
    {6,  {1,0,0,0,0,1,1}},
    {5,  {1,0,1,0,1,1}},

    //0x9
    {4,  {1,0,0,1,1}},
    {3,  {1,0,1,1}},
    {2,  {1,1,1}},
    {1,  {1,1}}
    };
        

    if(crc_size != 0){
    _crc_size = crc_size;
    _crc_array.resize(crc_size+1);
    std::vector<uint8_t> crc = crc_map[crc_size];
    _crc_array = crc;
    }else if(crc_size==0){
        _crc_size = 0;
        _crc_array.resize(0);
    }
    else{
        std::cout << "CRC size not supported" << std::endl;
        exit(0);
    }
}


void PACcode::create_bit_rev_order() {
    for (uint16_t i = 0; i < _block_length; ++i) {
        uint16_t to_be_reversed = i;
        _bit_rev_order.at(i) = (uint16_t) ((to_be_reversed & 1) << (_n - 1));
        for (uint8_t j = (uint8_t) (_n - 1); j; --j) {
            to_be_reversed >>= 1;
            _bit_rev_order.at(i) += (to_be_reversed & 1) << (j - 1);
        }
    }
}

//PAC code rate profile
std::vector<uint8_t> PACcode::rate_profile(std::vector<uint8_t> info_bits){
    std::vector<uint8_t> profiled_bits(_block_length);

    int count = 0;
    for(int i=0;i<_block_length;i++){
        if(_frozen_bits.at(i) == 0){
            profiled_bits.at(i) = info_bits.at(count);
            count++;
        }
        else{
            profiled_bits.at(i) = 0;
        }
    }

    return profiled_bits;
}

std::vector<uint8_t> PACcode::convolutional_encoder(std::vector<uint8_t> profiled_bits){

    std::vector<uint8_t> convoluted_bits(_block_length,0);

    //convolutional operation between profiled bits and convolutional polynomial
    for(int i = 0;i<_block_length;i++){
        for(int j=0;j<_m;j++){
            if(i-j>=0){
                convoluted_bits.at(i) = convoluted_bits.at(i) ^ (profiled_bits.at(i-j) & convolutional_polynomial.at(j));
            }
        }        
    }

    return convoluted_bits;
}


std::vector<uint8_t> PACcode::polar_encoder(std::vector<uint8_t> info_bits) {

    std::vector<uint8_t> coded_bits(_block_length);

    for (uint8_t iteration = 0; iteration < _n; ++iteration) {
        uint16_t  increment = (uint16_t) (1 << iteration);
        for (uint16_t j = 0; j < increment; j +=  1) {
            for (uint16_t i = 0; i < _block_length; i += 2 * increment) {
                info_bits.at(i + j) = (uint8_t)((info_bits.at(i + j) + info_bits.at(i + j + increment)) % 2);
            }
        }
    }

    for (uint16_t i = 0; i < _block_length; ++i) {
        coded_bits.at(i) = info_bits.at(_bit_rev_order.at(i));
    }

    return coded_bits;

}

std::vector<uint8_t> PACcode::decode_scascl_llr(std::vector<double> llr, uint16_t list_size) {

    _list_size = list_size;

    _llr_based_computation = true;

    initializeDataStructures();

    uint16_t  l = assignInitialPath();

    double * llr_0 = getArrayPointer_LLR(0, l);

    for (uint16_t beta = 0; beta < _block_length; ++beta ) {
        llr_0[beta] = llr.at(beta);
    }

    return decode_scascl();

}


std::vector<uint8_t> PACcode::decode_scascl() {
    
    uint16_t unfrozenbit_count = 0;
    bool early_termination = true;
    uint16_t info_crc_index;
    int phi2=0;
    
    _current_segment = 0;
    _segmented_block_first_index = 0;
    _current_list_size = _segmented_list_size.at(0);

    //initialize _shift_register
    for(uint16_t l=0;l<_max_list_size;l++){
        for(uint16_t i=0;i<_m;i++){
            _shift_register.at(l).at(i) = 0;
        }
    }
    
    //main loop
    for (uint16_t phi = 0; phi < _block_length; ++phi ){

        recursivelyCalcLLR(_n, phi);

        if (_frozen_bits.at(phi) == 1){
            continuePaths_FrozenBit(phi);
        }
        else{
            continuePaths_UnfrozenBit(phi);

            unfrozenbit_count++;
        }

        if ((phi%2) == 1){
            recursivelyUpdateC(_n, phi);
        }
    }
    

    std::vector<uint8_t> decoded_info_bits(_info_length,0); 

    uint16_t l = findMostProbablePath((bool) true); 

    uint8_t * c_0 = _arrayPointer_Info.at(l);

    std::vector<uint8_t> info_crc_bits(_info_length + _crc_mode);

    info_crc_index = 0;

    for (uint16_t beta = 0; beta < _block_length; ++beta ){

        if(_frozen_bits.at(beta) == 0){
            info_crc_bits.at(info_crc_index) = c_0[beta];
            info_crc_index++;
        }

        if(info_crc_index == _info_length + _crc_mode){
            break;
        }
    }

    uint16_t info_index=0;
    info_crc_index = 0;

    for(int s=0;s<_S;s++){
        for(int i=0;i<_segmented_bits_length.at(s) - _segmented_crc_length.at(s);i++){

            decoded_info_bits.at(info_index) = info_crc_bits.at(info_crc_index);
            info_index++;
            info_crc_index++;
        }
        info_crc_index += _segmented_crc_length.at(s);
    }

    for (uint16_t s = 0; s < _max_list_size; ++s) {
        delete[] _arrayPointer_Info.at(s);
        for (uint16_t lambda = 0; lambda < _n + 1; ++lambda) {

            delete[] _arrayPointer_LLR.at(lambda).at(s);
            delete[] _arrayPointer_C.at(lambda).at(s);
        }
    }

    return decoded_info_bits;

}



void PACcode::initializeDataStructures() {

    while (_inactivePathIndices.size()) {
        _inactivePathIndices.pop();
    };
    _activePath.resize(_max_list_size);

    _pathMetric_LLR.resize(_max_list_size);
    _arrayPointer_LLR.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _arrayPointer_LLR.at(i).resize(_max_list_size);

    _arrayPointer_C.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _arrayPointer_C.at(i).resize(_max_list_size);

    _arrayPointer_Info.resize(_max_list_size);

    _pathIndexToArrayIndex.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _pathIndexToArrayIndex.at(i).resize(_max_list_size);

    _inactiveArrayIndices.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i) {
        while (_inactiveArrayIndices.at(i).size()) {
            _inactiveArrayIndices.at(i).pop();
        };
    }

    _arrayReferenceCount.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _arrayReferenceCount.at(i).resize(_max_list_size);

    for (uint16_t s = 0; s < _max_list_size; ++s) {
        _arrayPointer_Info.at(s) = new uint8_t[_block_length]();
        for (uint16_t lambda = 0; lambda < _n + 1; ++lambda) {

            _arrayPointer_LLR.at(lambda).at(s) = new double[(1 << (_n - lambda))]();

            _arrayPointer_C.at(lambda).at(s) = new uint8_t[2 * (1 << (_n - lambda))]();
            _arrayReferenceCount.at(lambda).at(s) = 0;
            _inactiveArrayIndices.at(lambda).push(s);
        }
    }

    for (uint16_t l = 0; l < _max_list_size; ++l) {
        _activePath.at(l) = 0;
        _inactivePathIndices.push(l);
        _pathMetric_LLR.at(l) = 0;
    }
}

uint16_t PACcode::assignInitialPath() {

    uint16_t  l = _inactivePathIndices.top();
    _inactivePathIndices.pop();
    _activePath.at(l) = 1;
    // Associate arrays with path index
    for (uint16_t lambda = 0; lambda < _n + 1; ++lambda) {
        uint16_t  s = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();
        _pathIndexToArrayIndex.at(lambda).at(l) = s;
        _arrayReferenceCount.at(lambda).at(s) = 1;
    }
    return l;
}

uint16_t PACcode::clonePath(uint16_t l) {

    uint16_t l_p = _inactivePathIndices.top();
    _inactivePathIndices.pop();
    _activePath.at(l_p) = 1;

    _pathMetric_LLR.at(l_p) = _pathMetric_LLR.at(l);

    for (uint16_t lambda = 0; lambda < _n + 1; ++lambda ) {
        uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
        _pathIndexToArrayIndex.at(lambda).at(l_p) = s;
        _arrayReferenceCount.at(lambda).at(s)++;
    }
    return l_p;
}

void PACcode::killPath(uint16_t l) {

    _activePath.at(l) = 0;
    _inactivePathIndices.push(l);

    _pathMetric_LLR.at(l) = 0;

    for (uint16_t lambda = 0; lambda < _n + 1; ++lambda ) {
        uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
        _arrayReferenceCount.at(lambda).at(s)--;

        if (_arrayReferenceCount.at(lambda).at(s) == 0 ) {
            _inactiveArrayIndices.at(lambda).push(s);
        }
    }
}

double * PACcode::getArrayPointer_P(uint16_t lambda, uint16_t  l) {

    uint16_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint16_t s_p;

    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {

        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        std::copy(_arrayPointer_P.at(lambda).at(s), _arrayPointer_P.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_P.at(lambda).at(s_p));
        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;
    }
    return _arrayPointer_P.at(lambda).at(s_p);
}

double * PACcode::getArrayPointer_LLR(uint16_t lambda, uint16_t  l) {

    uint16_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint16_t s_p;

    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {
        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));
        std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) +  (1 << (_n - lambda)),  _arrayPointer_LLR.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;
    }

    return _arrayPointer_LLR.at(lambda).at(s_p);
}


uint8_t * PACcode::getArrayPointer_C(uint16_t lambda, uint16_t  l) {

    uint16_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint16_t s_p;

    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {
        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) +  (1 << (_n - lambda)),  _arrayPointer_LLR.at(lambda).at(s_p));
        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;

    }
    return _arrayPointer_C.at(lambda).at(s_p);
}

void PACcode::recursivelyCalcP(uint16_t lambda, uint16_t phi) {
    if ( lambda == 0 )
        return;

    uint16_t psi = phi >> 1;

    if ( (phi % 2) == 0)

        recursivelyCalcP(lambda -1, psi);

    double sigma = 0.0f;

    for (uint16_t l = 0; l < _max_list_size; ++l) {
        if (_activePath.at(l) == 0){
            continue;
        }

        double * p_lambda = getArrayPointer_P(lambda, l);
        double * p_lambda_1 = getArrayPointer_P(lambda - 1, l);
        
        uint8_t * c_lambda = getArrayPointer_C(lambda, l);
        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {

            if ( (phi %2) == 0 ){

                p_lambda[2 * beta] = 0.5f * ( p_lambda_1[2*(2*beta)]*p_lambda_1[2*(2*beta+1)]
                                              + p_lambda_1[2*(2*beta) + 1]*p_lambda_1[2*(2*beta+1) + 1]);

                p_lambda[2 * beta + 1] = 0.5f * ( p_lambda_1[2*(2*beta) +1]*p_lambda_1[2*(2*beta+1)]
                                                  + p_lambda_1[2*(2*beta)]*p_lambda_1[2*(2*beta+1) + 1]);
            }
            else {

                uint8_t  u_p = c_lambda[2*beta];

                p_lambda[2 * beta] = 0.5f * p_lambda_1[2*(2*beta) + (u_p % 2)] *   p_lambda_1[2*(2*beta + 1)];

                p_lambda[2 * beta + 1] = 0.5f * p_lambda_1[2*(2*beta) + ((u_p+1) % 2)] *   p_lambda_1[2*(2*beta + 1) + 1];
            }

            sigma = std::max(sigma,  p_lambda[2 * beta]);
            sigma = std::max(sigma,  p_lambda[2 * beta + 1]);


        }
    }

    for (uint16_t l = 0; l < _max_list_size; ++l) {

        if (sigma == 0) // Typically happens because of undeflow
            break;

        if (_activePath.at(l) == 0){
            continue;
        }

        double *p_lambda = getArrayPointer_P(lambda, l);
        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {

            p_lambda[2 * beta] = p_lambda[2 * beta] / sigma;
            p_lambda[2 * beta + 1] = p_lambda[2 * beta + 1] / sigma;
        }
    }
}

void PACcode::recursivelyCalcLLR(uint16_t lambda, uint16_t phi) {

    if ( lambda == 0 )
        return;

    uint16_t psi = phi >> 1;

    if ( (phi % 2) == 0){
        recursivelyCalcLLR(lambda -1, psi);
    }

    for (uint16_t l = 0; l < _max_list_size; ++l) {

        if (_activePath.at(l) == 0){
            continue;
        }

        double * llr_lambda = getArrayPointer_LLR(lambda, l);
        double * llr_lambda_1 = getArrayPointer_LLR(lambda - 1, l);

        uint8_t * c_lambda = getArrayPointer_C(lambda, l);
        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {

            if ( (phi %2) == 0 ){

                if (40 > std::max(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]))){

                    llr_lambda[beta] = std::log ( (exp(llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1]) + 1) /
                                                  (exp(llr_lambda_1[2*beta]) + exp(llr_lambda_1[2*beta+1])));
                }
                else {

                    llr_lambda[beta] = (double)  ((llr_lambda_1[2 * beta] < 0) ? -1 : (llr_lambda_1[2 * beta] > 0)) *
                                       ((llr_lambda_1[2 * beta + 1] < 0) ? -1 : (llr_lambda_1[2 * beta + 1] > 0)) *
                                       std::min( std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]));
                }
            }
            else {
                
                uint8_t  u_p = c_lambda[2*beta];

                llr_lambda[beta] = (1 - 2 * u_p) * llr_lambda_1[2*beta] + llr_lambda_1[2*beta + 1];
            }
        }
    }
}

void PACcode::recursivelyUpdateC(uint16_t lambda, uint16_t phi) {
    uint16_t psi = phi >> 1;
    
    for (uint16_t l = 0; l < _max_list_size; ++l) {
        if (_activePath.at(l) == 0){
            continue;
        }

        uint8_t *c_lambda = getArrayPointer_C(lambda, l);
        uint8_t *c_lambda_1 = getArrayPointer_C(lambda - 1, l);

        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {

            c_lambda_1[2 * (2 * beta) + (psi % 2)] = (uint8_t) ((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);

            c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
        }
    }

    if ( (psi % 2) == 1){
        recursivelyUpdateC((uint16_t) (lambda - 1), psi);
    }
}


void PACcode::continuePaths_FrozenBit(uint16_t phi) {

    int u_hat;

    std::vector<uint8_t>  v(_m);
    
    for (uint16_t l = 0; l < _max_list_size; ++ l) {
        if (_activePath.at(l) == 0){
            continue;
        }


        //left shift _shift_register, with the rightmost bit set to 0
        for(int j=0;j<_m-1;j++){
            _shift_register.at(l).at(j) = _shift_register.at(l).at(j+1);
        }
        _shift_register.at(l).at(_m-1) = 0;

        uint8_t  * c_m = getArrayPointer_C(_n, l);

        //v <= _shift_register
        for(int j=0;j<_m;j++){
            v.at(j) = _shift_register.at(l).at(j);
        }

        u_hat = 0;
        for(int j=0;j<_m;j++){
            if(phi-j>=0){
                u_hat = u_hat ^ (convolutional_polynomial.at(j) & v.at(_m-1-j));
            }
        }

        // printf("%d",u_hat);

        c_m[(phi % 2)] = u_hat; 
        // c_m[(phi % 2)] = 0; 

        double *llr_p = getArrayPointer_LLR(_n, l); 

        _pathMetric_LLR.at(l) += log(1 + exp(-(1.0-2.0*u_hat)*llr_p[0]));

        _arrayPointer_Info.at(l)[phi] = 0;
    }
}



void PACcode::continuePaths_UnfrozenBit(uint16_t phi) {

    std::vector<double>  probForks((unsigned long) (2 * _max_list_size));
    std::vector<double> probabilities;
    std::vector<uint8_t>  contForks((unsigned long) (2 * _max_list_size));

    int u_hat, u_hat_p;
    std::vector<uint8_t> v(_m);

    uint16_t  i = 0;
    for (unsigned l = 0; l < _max_list_size; ++l) {

        if (_activePath.at(l) == 0) {
            probForks.at(2 * l) = NAN;
            probForks.at(2 * l + 1) = NAN;
        }
        else {
            double *llr_p = getArrayPointer_LLR(_n, l);

            probForks.at(2 * l) =  - (_pathMetric_LLR.at(l) + log(1 + exp(-llr_p[0])));

            probForks.at(2 * l + 1) = -  (_pathMetric_LLR.at(l) + log(1 + exp(llr_p[0])));

            probabilities.push_back(probForks.at(2 * l));
            probabilities.push_back(probForks.at(2 * l +1));

            i++;
        }
    }

    uint16_t  rho = _current_list_size;

    if ( (2*i) < _current_list_size){
        rho = (uint16_t) 2 * i;
    }


    for (uint8_t l = 0; l < 2 * _max_list_size; ++l) {
        contForks.at(l) = 0;
    }

    std::sort(probabilities.begin(), probabilities.end(), std::greater<double>());

    double threshold = probabilities.at((unsigned long) (rho - 1));

    uint16_t num_paths_continued = 0;

    for (uint8_t l = 0; l < 2 * _max_list_size; ++l) {

        if (probForks.at(l) > threshold) {
            contForks.at(l) = 1;
            num_paths_continued++;
        }
        if (num_paths_continued == rho) {
            break;
        }
    }

    if  ( num_paths_continued < rho ) {
        for (uint8_t l = 0; l < 2 * _max_list_size; ++l) {

            if (probForks.at(l) == threshold) {
                contForks.at(l) = 1;
                num_paths_continued++;
            }

            if (num_paths_continued == rho) {
                break;
            }
        }
    }

    for (unsigned l = 0; l < _max_list_size; ++l) {

        if (_activePath.at(l) == 0){
            continue;
        }

        if ( contForks.at(2 * l)== 0 && contForks.at(2 * l + 1) == 0 ){
            killPath(l);
        }
    }

    for (unsigned l = 0; l < _max_list_size; ++l) {

        if ( contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0 ){
            continue;
        }

        uint8_t * c_m = getArrayPointer_C(_n, l);

        //left-shift _shift_register, with the rightmost bit set to 0
        for(int j=0;j<_m-1;j++){
            _shift_register.at(l).at(j) = _shift_register.at(l).at(j+1);
        }
        _shift_register.at(l).at(_m-1) = 0;

        //v <= _shift_register
        for(int j=0;j<_m;j++){
            v.at(j) = _shift_register.at(l).at(j);
        }

        if ( contForks.at(2 * l) == 1 && contForks.at(2 * l + 1) == 1 ) {
            
            u_hat = 0;
            for(int j=0;j<_m;j++){
                if(phi-j>=0){
                    u_hat = u_hat ^ (convolutional_polynomial.at(j) & v.at(_m-1-j));
                }
            }

            c_m[(phi%2)] = u_hat;

            uint16_t l_p = clonePath(l);
            c_m = getArrayPointer_C(_n, l_p);

            //_shift_register[l_p] <= _shift_register[l]
            for(int j=0;j<_m;j++){
                _shift_register.at(l_p).at(j) = _shift_register.at(l).at(j);
            }

            //flip the rightmost bit of _shift_register
            _shift_register.at(l_p).at(_m-1) = _shift_register.at(l_p).at(_m-1) ^ 1;

            //v<=_shift_register[l_p]
            for(int j=0;j<_m;j++){
                v.at(j) = _shift_register.at(l_p).at(j);
            }

            u_hat_p = 0;
            for(int j=0;j<_m;j++){
                if(phi-j>=0){
                    u_hat_p = u_hat_p ^ (convolutional_polynomial.at(j) & v.at(_m-1-j));
                }
            }

            c_m[(phi%2)] = u_hat_p;

            std::copy(_arrayPointer_Info.at(l), _arrayPointer_Info.at(l) +  phi,  _arrayPointer_Info.at(l_p));

            _arrayPointer_Info.at(l)[phi] = _shift_register.at(l).at(_m-1);
            _arrayPointer_Info.at(l_p)[phi] = _shift_register.at(l_p).at(_m-1);

            double *llr_p = getArrayPointer_LLR(_n, l);
            _pathMetric_LLR.at(l) += log(1 + exp(-(1.0-2.0*u_hat)*llr_p[0]));

            llr_p = getArrayPointer_LLR(_n, l_p);
            _pathMetric_LLR.at(l_p) += log(1 + exp(-(1.0-2.0*u_hat_p)*llr_p[0]));
            
        }
        else {
            if ( contForks.at(2 * l) == 1) {
                u_hat = 0;
                for(int j=0;j<_m;j++){
                    if(phi-j>=0){
                        u_hat = u_hat ^ (convolutional_polynomial.at(j) & v.at(_m-1-j));
                    }
                }

                if(u_hat==1){
                    //flip the rightmost bit of _shift_register
                    _shift_register.at(l).at(_m-1) = _shift_register.at(l).at(_m-1) ^ 1;
                }

                c_m[(phi%2)] = 0;

                _arrayPointer_Info.at(l)[phi] = _shift_register.at(l).at(_m-1);

                double *llr_p = getArrayPointer_LLR(_n, l);
                _pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
                
            }
            else {
                u_hat = 0;
                for(int j=0;j<_m;j++){
                    if(phi-j>=0){
                        u_hat = u_hat ^ (convolutional_polynomial.at(j) & v.at(_m-1-j));
                    }
                }

                if(u_hat==0){
                    //flip the rightmost bit of _shift_register
                    _shift_register.at(l).at(_m-1) = _shift_register.at(l).at(_m-1) ^ 1;
                }

                c_m[(phi%2)] = 1;

                _arrayPointer_Info.at(l)[phi] = _shift_register.at(l).at(_m-1);

                double *llr_p = getArrayPointer_LLR(_n, l);
                _pathMetric_LLR.at(l) += log(1 + exp(llr_p[0]));
                
            }
        }
    }

}

uint16_t PACcode::findMostProbablePath(bool check_crc) {

    uint16_t  l_p = 0;
    double p_p1 = 0;
    double p_llr = std::numeric_limits<double>::max();
    bool path_with_crc_pass = false;


    for (uint16_t l = 0; l < _max_list_size; ++l) {


        if (_activePath.at(l) == 0){
            continue;
        }

        path_with_crc_pass = true;

        if (_pathMetric_LLR.at(l) < p_llr ) {

            p_llr = _pathMetric_LLR.at(l);
            l_p  = l;
        }
        
    }
    
    if ( path_with_crc_pass)
        return l_p;
    else
        return findMostProbablePath(false);
}

double PACcode::calcXi(double x){
  double k0, y;
  double a, b, c;
  double X1, X2, X3, Y1, Y2;
  double a0, a1, a2;
  
  X1 = 0.2; X2 = 0.7; X3 = 10.0;
  Y1 = 0.125 * pow( X1, 2.0 ) - 0.5 * X1 - 0.125 * pow( X1, 3.0 );
  a = - 0.4527; b = 0.0218; c = 0.86;
  Y2 = a * pow( X2, c ) + b;
  
  if( x <= X1 )
     y = 0.125 * pow( x, 2.0 ) - 0.5 * x - 0.125 * pow( x, 3.0 );
  else if( x <= X2 ){
    a0 = -0.002706; a1 = -0.476711; a2 = 0.051200;

    y = a0 + a1 * x + a2 * pow( x, 2.0 );
  }  else if( x < X3 ){
    // Chung
    a = - 0.4527; b = 0.0218; c = 0.86;
    y = a * pow( x, c ) + b;
  }  else  {
    k0 = 8.554;
    y =  - 0.25 * x - 0.5 * log( x ) + 0.5 * log( M_PI ) + log( 1.0 - M_PI * M_PI / (4.0 * x )  + k0 / ( x * x ) );
  }

  return y;
}

double PACcode::calcXiInverse(double y){

  double k0, x,eps;
  double a, b, c;
  double X1, X2, X3, Y1, Y2, Y3;
  double a0, a1, a2;
  X1 = 0.2; X2 = 0.7; X3 = 10.0;

  Y1 = 0.125 * powl( X1, 2.0 ) - 0.5 * X1 - 0.125 * powl( X1, 3.0 );

  a0 = -0.002706; a1 = -0.476711; a2 = 0.051200;
  Y2 = a0 + a1 * X2 + a2 * powl( X2, 2.0 );

  // Chung
  a = - 0.4527; b = 0.0218; c = 0.86;
  Y3 = a * powl( X3, c ) + b;

  if( y >= Y1 ) 
    {
      x = - 2.0 * y + powl( y, 2.0 ) + powl( y, 3.0 );
    } 
  else if( y>= Y2 ){
    x = ( - a1 - sqrt( a1 * a1 - 4.0 * a2 * ( a0 - y ) ) ) / (2.0 * a2 );
  } else if( y >= Y3 ){
    a = - 0.4527; b = 0.0218; c = 0.86;
    x = powl( ( y - b ) / a, 1.0 / c );
  } else {
    eps = 1.0e-12;
    if( y < -2000.0 )
      eps = 1.0e-10;

    x = bisection( y, 0.0, -4.0 * y, eps );
    if( x < -1.0 ) {
      printf("retry\n");
      eps = 1.0e-9;
      x = bisection( y, 0.0, -4.0 * y, eps );      
    }
    if( x < -1.0 ) {
      printf("retry again\n");
      eps = 1.0e-8;
      x = bisection( y, 0.0, -4.0 * y, eps );      
    }
    if( x < -1.0 ) {
      printf("Give up\n");
      exit(0);
    }

  }

  return x;
}

double PACcode::bisection(double y, double a, double b, double eps ){
    double s, diff;
    int i=0;

    while (!(fabs(a-b)<eps)){
        if( i > 10000 ) {
    diff = fabs( a- b );
    printf("failed to converge., y = %f, a = %f, b = %f, diff = %f(%e)\n", y, a, b, diff, diff );
    return - 100.0;
        }
        i++;

        s = (a+b)/2.0;
        if( ff(s,y) * ff(a,y)<0.0) b=s;
        else a = s;

    };
    return s;    
}

double PACcode::ff(double x, double y)
{
  double f;
  double k0;


  k0 = 8.554;
  f =  - 0.25 * x - 0.5 * log( x ) + 0.5 * log( M_PI ) + log( 1.0 - M_PI * M_PI / (4.0 * x )  + k0 / ( x * x ) ) - y;

  return f;
}

double PACcode::calcGamma(double x){
  double w, y, z;
  double X1, X2, X3, Y1, Y2;
  X1 = 0.2; X2 = 0.7; X3 = 10.0;
  double c0, c1, c2, c3, c4, c5;


  if( x <= X1 ){
    z = 0.5 * pow( x, 2.0 ) - 0.5 * pow( x, 3.0 ) + 2.0 / 3.0 * pow( x, 4.0 );
  }
  else{
    y = calcXi( x );
    z = calcXiInverse( y + log( 2.0 - exp( y ) ) );
  }

  return z;
}

double PACcode::CALC_CAP(double xi){
    double A;
    double B;
    double U;
    double gamma;

    double alpha = 1.16125;
    double threshold1 = 0.04;
    double threshold2 = 1.0;
    double threshold3 = 10.0;
    double big_xi = -11.3143;

    double C1 = 0.055523;
    double C2 = 0.721452;

    double H21 = 1.396634;
    double H22 = 0.872764;
    double H23 = 1.148562;
    double H31 = 1.266967;
    double H32 = 0.938175;
    double H33 = 0.986830;

    double ln2 = log(2);
    if(xi < big_xi){
        B = ln2 + 2.0 * log(ln2) + 2.0 * log(alpha) - 2.0 * xi;
        
        return log(B + (1.0/B - 1.0) * log(B)) - ln2;
    }

    gamma = exp(xi);

    if(gamma > threshold3){
        return log(ln2) + log(alpha) - gamma - 0.5 * xi;
    }
    else if(gamma < threshold1){
        U = 1.0 - (gamma - gamma * gamma + 4.0 / 3.0 * pow(gamma,3.0))/ln2;
    }
    else if(gamma < threshold2){
        U = 1.0 - pow(1.0 - exp(- H21 * pow(gamma,H22)),H23);
    }
    else{
        U = 1.0 - pow(1.0 - exp(- H31 * pow(gamma,H32)),H33);
    }

    if(U < C1){
        
        
        A = pow(-5.0 + 24.0 * ln2 * U + 2.0 * sqrt(13.0 + 12.0 * ln2 * U * (12.0 * ln2 * U - 5.0)), 1.0/3.0);
        
        
        return log(1.0 - 3.0 / A + A) - 2.0 * ln2;
    }
    else if(U < C2){

        return ( log(-log(1.0 - pow(U , 1.0 / H23))) - log(H21)) / H22;
    }
    else{
        return ( log(-log(1.0 - pow(U , 1.0 / H33))) - log(H31)) / H32;
    }
}

void PACcode::initialize_frozen_bits() {
    std::vector<double> channel_vec(_block_length);

    if(_construction_mode == 0){
        printf("construction mode is BTC\n");

        channel_vec.at(0) = exp(-pow(10.0,0.1 * _design_snr));
        for(int j = 0; j < _n; j ++ ){
        int u = ( 1 << (j+1) );

            for(int t = 0; t < (u/2); t++ ){
                double T = channel_vec.at(t);
                channel_vec.at(t) = 2.0 * T - T * T;
                channel_vec.at(t+u/2) = T * T;
            }
        }

        _channel_order_descending.resize(_block_length);
        std::size_t n_t(0);
        std::generate(std::begin(_channel_order_descending), std::end(_channel_order_descending), [&]{ return n_t++; });
        std::sort(  std::begin(_channel_order_descending),
                    std::end(_channel_order_descending),
                    [&](int i1, int i2) { return channel_vec[_bit_rev_order.at(i1)] < channel_vec[_bit_rev_order.at(i2)]; } );
    }
    else if(_construction_mode == 1){

        printf("construction mode is IGA\n");
        channel_vec.at(0) = 4.0 * pow(10.0,0.1 * _design_snr);

        uint16_t J;
        double u;

        for(int i=1;i<=_n;i++){
            uint16_t J = (1 << i);

            for(int j=0;j<=(int)(J/2)-1;j++){
                u = channel_vec.at(j);
                channel_vec.at(j) = calcGamma(u);
                channel_vec.at(j + (int)J/2) = 2.0 * u;
            }
        }


        _channel_order_descending.resize(_block_length);
        std::size_t n_t(0);
        std::generate(std::begin(_channel_order_descending), std::end(_channel_order_descending), [&]{ return n_t++; });
        std::sort(  std::begin(_channel_order_descending),
                    std::end(_channel_order_descending),
                    [&](int i1, int i2) { return channel_vec[_bit_rev_order.at(i1)] > channel_vec[_bit_rev_order.at(i2)]; } );
     
    }
    else if(_construction_mode == 2){

        printf("construction mode is RCA\n");
        channel_vec.at(0) = log(pow(10.0,0.1 * _design_snr));

        
        double xi0;
        double cap0;
        int J;

        double ln2 = log(2);

        for(int i=1;i<=_n;i++){
            J = (1 << i);
            
            for(int j=0; j<= (int)J/2 -1 ;j++){

                xi0 = channel_vec.at(j);
                cap0 = CALC_CAP(xi0);

                channel_vec.at(j) = CALC_CAP(cap0 + ln2);
                channel_vec.at(j + (int)J/2) = xi0 + ln2;

            }
        }

        _channel_order_descending.resize(_block_length);
        std::size_t n_t(0);
        std::generate(std::begin(_channel_order_descending), std::end(_channel_order_descending), [&]{ return n_t++; });
        std::sort(  std::begin(_channel_order_descending),
                    std::end(_channel_order_descending),
                    [&](int i1, int i2) { return channel_vec[_bit_rev_order.at(i1)] > channel_vec[_bit_rev_order.at(i2)]; } );
    }
    else if(_construction_mode == 3){

        printf("construction mode is RM-profile\n");

        //wt(i), which is the number of 1s in the binary representation of i
        //information bits are chosen as k bits with the largest wt(i)
        //frozen bits are chosen as n-k bits with the smallest wt(i)

        std::vector<uint16_t> wt(_block_length);
        std::vector<uint16_t> wt_sorted(_block_length);

        for(int i=0;i<_block_length;i++){
            wt.at(i) = 0;
            // wt_sorted.at(i) = 0;
        }

        for(int i=0;i<_block_length;i++){
            uint16_t temp = i;
            while(temp > 0){
                wt.at(i) += temp % 2;
                temp = temp / 2;
            }
        }

        _channel_order_descending.resize(_block_length);
        std::size_t n_t(0);
        std::generate(std::begin(_channel_order_descending), std::end(_channel_order_descending), [&]{ return n_t++; });
        std::sort(  std::begin(_channel_order_descending),
                    std::end(_channel_order_descending),
                    [&](int i1, int i2) { return wt[_bit_rev_order.at(i1)] > wt[_bit_rev_order.at(i2)]; } );

    }
    else{
        exit(0);
    }
 
    uint16_t  effective_info_length = _info_length + _crc_mode;

    for (uint16_t i = 0; i < effective_info_length; ++i) {
        _frozen_bits.at(_channel_order_descending.at(i)) = 0;
    }
    for (uint16_t i = effective_info_length; i < _block_length; ++i) {
        _frozen_bits.at(_channel_order_descending.at((i))) = 1;
    }
}


uint8_t PACcode::segment_crc_check(std::vector<uint8_t> info_bit_padded) {

    std::vector<uint8_t> crc_pool(_segmented_bits_length.at(_current_segment));
    bool early_termination = false;
    uint8_t pass_check = 1;

    crc_type(_segmented_crc_length.at(_current_segment));
    uint16_t info_crc_index = 0;

    for(int i = 0; i< _segmented_bits_length.at(_current_segment); ++i){
        crc_pool.at(i) = info_bit_padded.at(i);
    }

    for(int i=0;i<_segmented_bits_length.at(_current_segment) - _segmented_crc_length.at(_current_segment);i++){

        if(crc_pool.at(i) != 0){
            for(int j=0;j < _segmented_crc_length.at(_current_segment)+1;j++){

                crc_pool.at(i+j) = crc_pool.at(i+j)^ _crc_array.at(j);
            }
        }
    }

    for (uint16_t i = _segmented_bits_length.at(_current_segment) - _segmented_crc_length.at(_current_segment); i < _segmented_bits_length.at(_current_segment); ++i) {

        if(crc_pool.at(i) != 0){
            early_termination = false;
            pass_check = 0;
            break;
        }
    }

    return pass_check;
}




