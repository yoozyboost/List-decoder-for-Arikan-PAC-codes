#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <stack>
#include <set>
#include <sstream>
#include <string>
#include <iterator>
#include <fstream>
#include "PACcode.h"

using namespace std;

#define N 64
// #define r 
// #define k (512 - r)
#define l 2
#define S 1//number of segments
// #define CODE_RATE (double)k/N
#define n (int)log2(N)

#define PI M_PI

#define SNR_MIN -3				
#define SNR_MAX 3
#define EbN0_MIN 2.0
#define EbN0_MAX 2.0
#define SNR_INTERVAL 1
#define EbN0_INTERVAL 0.5
#define TRIAL_MAX 1000000
#define TRIAL_MIN 1000
#define TRIAL_INTERVAL 1000		
#define ERROR_MAX 1000



int main(int argc, char *argv[]){

    long int trial;
	double snr;//signal noise ratio	
    double ebn0;				    
    double BER; //bit error rate
    double BLER; //block error rate
    std::vector<double> segmented_BLER(S);
    std::vector<double> segmented_after_BLER(S);

    double EbN0dB;
    double EbN0;
    double EsN0;

    double Eb;
    double Es;
    double N0;

    double sigma_w; //variance of AWGN

    int error_count;
    int block_error_count;

    int bound_after_priority_block_error_count;
    int after_priority_block_error_count;
    int priority_block_error_count;
    int regular_block_error_count;

    std::vector<int> segmented_block_error_count(S);
    std::vector<int> segmented_after_block_error_count(S);

    

    int construction_mode;
    double design_snr;
    int crc_mode;

    double average_list_size;

    srand((unsigned) time(NULL));
    
    int k;
    int list;
    int k0;
    int r;

    design_snr = -1.0;

    construction_mode = 3;
    
    r = 0;
    k = 32;

    double CODE_RATE = (double)k/N;
    
    int flag;
    int error_count_to_print;    

    std::vector<uint16_t> segmented_list_size(S);
    for (int i = 0; i < S; i++)
    {
        segmented_list_size.at(i) = l;
    }
    
    std::vector<uint16_t> segmented_bits_length(S);
    std::vector<uint16_t> segmented_crc_length(S);

    segmented_bits_length.at(0) = k;
    segmented_crc_length.at(0) = 0;
    int total_r = r;


    int list_size = *std::max_element(segmented_list_size.begin(), segmented_list_size.end());

    // std::vector<std::string> sFile(2*S+1);
    // std::vector<std::string> sFile2(S);
    
    std::string sFile;
    std::stringstream ss1,ss2,ss3;


    for (auto it = segmented_bits_length.begin(); it != segmented_bits_length.end(); it++)    {
        if (it != segmented_bits_length.begin()) {
            ss1 << ",";
        }
        ss1 << *it;
    }

    for (auto it = segmented_list_size.begin(); it != segmented_list_size.end(); it++)    {
        if (it != segmented_list_size.begin()) {
            ss2 << ",";
        }
        ss2 << *it;
    }

    for (auto it = segmented_crc_length.begin(); it != segmented_crc_length.end(); it++)    {
        if (it != segmented_crc_length.begin()) {
            ss3 << ",";
        }
        ss3 << *it;
    }


    if(construction_mode == 0){
        // sFile1 = "result/SCASCL-BER-BPSK-AWGN-N=" + std::to_string(N) + "-R= " + std::to_string(CODE_RATE) + " -BTC-designSNR=" + std::to_string(design_snr) + "-L=" + std::to_string(list_size) + "-segment={" + ss1.str() +"}-crc={" + ss2.str() + "}.txt";
        // sFile2 = "result/SCASCL-BLER-BPSK-AWGN-N=" + std::to_string(N) + "-R= " + std::to_string(CODE_RATE) + " -BTC-designSNR=" + std::to_string(design_snr) + "-L=" + std::to_string(list_size) + "-segment={"+ ss1.str() +"}-crc={" + ss2.str() + "}.txt";
        // sFile3 = "result/SCASCL-averageListSize-BPSK-AWGN-N=" + std::to_string(N) + "-R= " + std::to_string(CODE_RATE) + " -BTC-designSNR=" + std::to_string(design_snr) + "-L=" + std::to_string(list_size) + "-segment={"+ ss1.str() +"}-crc={" + ss2.str() + "}.txt";

    }
    else if(construction_mode == 1){

        sFile = "result/PACcode-BLER-N=" + std::to_string(N) + "-R=" + std::to_string(CODE_RATE) + "-IGA-designSNR=" + std::to_string(design_snr) + "-segment={" + ss1.str() +"}-L={" + ss2.str() + "}-crc={" + ss3.str() + "}.txt";
        
    }
    else if(construction_mode == 2){

        // for(int i = 0; i < S; i++){
        //     sFile.at(i) = "result23/" + std::to_string(input_snr) + "@SCASCL-Segment" + std::to_string(i)+"_BLER-N=" + std::to_string(N) + "-R=" + std::to_string(CODE_RATE) + "-RCA-designSNR=" + std::to_string(design_snr) + "-segment={" + ss1.str() +"}-L={" + ss2.str() + "}-crc={" + ss3.str() + "}.txt";
        // }
        // for(int i=S;i<2*S;i++){
        //     sFile.at(i) = "result23/" + std::to_string(input_snr) + "@SCASCL-Segment" + std::to_string(i-S)+"_after_BLER-N=" + std::to_string(N) + "-R=" + std::to_string(CODE_RATE) + "-RCA-designSNR=" + std::to_string(design_snr) + "-segment={" + ss1.str() +"}-L={" + ss2.str() + "}-crc={" + ss3.str() + "}.txt";
        // }
        sFile = "result/PACcode-BLER-N=" + std::to_string(N) + "-R=" + std::to_string(CODE_RATE) + "-RCA-designSNR=" + std::to_string(design_snr) + "-segment={" + ss1.str() +"}-L={" + ss2.str() + "}-crc={" + ss3.str() + "}.txt";
    }

    std::ofstream fp;
    fp.open(sFile);
    
    PACcode paccode(N,k,construction_mode,design_snr,total_r,S,segmented_bits_length,segmented_crc_length,segmented_list_size);


    long int num_calc_count;
    long int priority_num_calc_count;

    std::vector<long double> partial_block_error_count(k,0);

    for( ebn0 = EbN0_MIN; ebn0 <= EbN0_MAX ; ebn0 += EbN0_INTERVAL ){
    // for( ebn0 = input_snr; ebn0 <= input_snr ; ebn0 += EbN0_INTERVAL ){

        num_calc_count = 0;
        priority_num_calc_count =0;

        sigma_w = sqrt(1) / sqrt( ( 
                            pow(10.0, ( (double)ebn0 / 10.0 ))));

        // Es = 1.0;
        // Eb = Es / (CODE_RATE); 
        // EbN0 = pow(10.0,0.1 * ebn0);
        // N0 = Eb / EbN0;
        // sigma_w = sqrt(N0);

		error_count = 0;
        block_error_count = 0;

        int check = 0;
        
        for( trial = 0 ; trial < TRIAL_MAX; trial++ ){


            std::vector<uint8_t> info_bits = paccode.dataGenerator(k);
            //show info_bits
            // for(int i=0;i<k;i++){
            //     printf("%d",info_bits.at(i));
            // }
            // printf("\n");
            std::vector<uint8_t> profiled_bits = paccode.rate_profile(info_bits);
            // show profiled_bits
            // for(int i=0;i<N;i++){
            //     printf("%d",profiled_bits.at(i));
            // }

            // printf("\n");
            std::vector<uint8_t> convoluted_bits = paccode.convolutional_encoder(profiled_bits);
            // //show convoluted_bits
            // for(int i=0;i<N;i++){
            //     printf("%d",convoluted_bits.at(i));
            // }

            // printf("\n");


        
            std::vector<uint8_t> code_bits = paccode.polar_encoder(convoluted_bits);
            // //show code_bits
            // for(int i=0;i<N;i++){
            //     printf("%d",code_bits.at(i));
            // }

            // printf("\n");

            // exit(0);
            std::vector<int8_t> sending_data = paccode.BPSKmodulation(code_bits);
            std::vector<double> received_signal = paccode.AWGNchannel(sending_data,sigma_w);

            std::vector<double> p0(N), p1(N);
            std::vector<double> llr(N);



            for (uint16_t i = 0; i < N; ++i) {
                llr.at(i) = 4 * received_signal.at(i) / (sigma_w * sigma_w);
            }


            std::vector<uint8_t> decoded_info_bits = paccode.decode_scascl_llr(llr, list_size);

            //check
            check += paccode.check(profiled_bits);

            // for(int i=0;i<k;i++){
            //     printf("%d",decoded_info_bits.at(i));
            // }
            // printf("\n");

            // exit(0);


            error_count += paccode.errorCount(info_bits,decoded_info_bits);
            block_error_count += paccode.blockErrorCount(info_bits,decoded_info_bits);            

            if( ( trial > TRIAL_MIN ) && ( block_error_count > ERROR_MAX ) ){
                trial++; break;
            }

            // show the progress
            if( trial % TRIAL_INTERVAL == 0 ){
            fprintf( stderr, " %f %8ld    (%d)\r", ebn0, trial, block_error_count);
            }

        }//end of one trial

    
        long all_bits = k * trial;
        long all_block = trial;
        BER = (double)error_count/all_bits;
        BLER = (double)block_error_count/all_block;

        double check_ratio = (double)check/all_block;

        // printf("check_ratio=%g\n",check_ratio);

        for(int i=0;i<S;i++){
            printf("EbN0=%g sigma=%g block error count = %d all block = %ld BLER=%g\n",ebn0,sigma_w,block_error_count,trial,BLER);
        }


        printf("\n\n");

        fp << ebn0 << " " << BLER << std::endl;

    }

    printf("finished!");

    return 0;
}