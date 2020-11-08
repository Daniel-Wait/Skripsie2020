#include <stdint.h>

#include "arm_math.h"
#include "arm_const_structs.h"

#define DOWN_SAMP 2
#define BIT_LEN 441 // SAMPLES PER SYMBOL ÷ DOWN_SAMP = 44.1kHz/20ms ÷ 2 = 441

#define L_LEN 2048
#define N_LEN 4096
#define B_BLOCKS 3
#define PACKETS 16 // 2048 SAMPLES ÷ 128 SAMPLES per PACKET

#include <Audio.h>
#include <Wire.h>
#include <SPI.h>
#include <SD.h>
#include <SerialFlash.h>

// GUItool: begin automatically generated code
AudioInputI2S            i2s2;
AudioMixer4              mixer1;
AudioOutputI2S           i2s1;
AudioRecordQueue         queue1;
AudioConnection          patchCord1(i2s2, 0, mixer1, 0);
AudioConnection          patchCord2(i2s2, 1, mixer1, 1);
AudioConnection          patchCord3(mixer1, 0, i2s1, 0);
AudioConnection          patchCord4(mixer1, 0, i2s1, 1);
AudioConnection          patchCord5(mixer1, queue1);
AudioControlSGTL5000     sgtl5000_1;
// GUItool: end automatically generated code

uint8_t bstart = 1;

int16_t fdata_head = 0;

int8_t detected_thresh = 0;
float corr_thresh = ((12+7)/2)*BIT_LEN/2;

const float f1 = 2200;
const float f0 = 1600;

float32_t mf_fft[B_BLOCKS][2*N_LEN] = {{0}};
float32_t rx_fft[B_BLOCKS][2*N_LEN] = {{0}};

const int ola_max = 3*N_LEN;
float32_t ola_buffer[ola_max] = {0};

int ola_head = N_LEN-L_LEN;
const int ola_overlap = N_LEN-L_LEN;

float32_t ola_peak = 0;
int ola_peak_index = 0;


/*SETUP*/
void setup() {
  // put your setup code here, to run once:
  float32_t mf_msg[B_BLOCKS][N_LEN] = {{0}};
  
  Serial.begin(300);
   
  AudioMemory(100);
  sgtl5000_1.enable();
  sgtl5000_1.inputSelect(AUDIO_INPUT_MIC);
  sgtl5000_1.micGain(0);
  sgtl5000_1.unmuteHeadphone();
  sgtl5000_1.audioPreProcessorEnable();
  sgtl5000_1.audioPostProcessorEnable();
  sgtl5000_1.volume(0.3);

  while (!Serial) {}

  const int bit_num = 12;
  
  int bitseq[bit_num] = {0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1};
  
  float space = f0*DOWN_SAMP;
  float mark = f1*DOWN_SAMP;
  generateTempate(mf_msg, bitseq, bit_num, 44100, BIT_LEN , space, mark);
  
  for(int p = 0; p < B_BLOCKS; p++)
  {
    f32tRealFFT(mf_fft[p], mf_msg[p], N_LEN);
  }
  
  queue1.begin();
}


/*LOOP*/
void loop() 
{
  if(bstart == 1)
  {
    if(detected_thresh < 1)
    {
      asm(" WFI");
      q15_t q_rx_msg[L_LEN*DOWN_SAMP] = {0};
       
      if (queue1.available() >= PACKETS*DOWN_SAMP)
      {  
        for (int k = 0; k < PACKETS*DOWN_SAMP; k++)
        {
          memcpy(q_rx_msg + 128*k, queue1.readBuffer(), 256);
          queue1.freeBuffer();
        }

        float32_t rx_msg[L_LEN*DOWN_SAMP] = {0};
        float32_t rx_down[N_LEN] = {0};
        
        arm_q15_to_float(q_rx_msg, rx_msg, L_LEN*DOWN_SAMP);
        downsample(rx_down, rx_msg, L_LEN);
            
        for(int k = 0; k < L_LEN; k++)  
        {
          Serial.println( rx_down[k] );
        }
        
        UPOLA(rx_down);
        check_detect();
      }   
    }
    else
    {
      /*RELEASE CONTROL*/
      
      queue1.clear();
      queue1.end();
      bstart = 0;
      
      /*PRINT*/
      for(int m = 0; m < ola_max; m++)
      {
        int indx = ola_head + m;
        if(indx >= ola_max)
        {
          indx -= ola_max;
        }            
        Serial.println( ola_buffer[indx] );
      }
      
    }
  }
}


/*DETECTION*/
void check_detect()
{
  if(ola_peak >= corr_thresh)
  {
    detected_thresh += 1;
    
    //int dist_LOW = BIT_LEN;
    float LOW_max_val = 0;

    for(int q = -10; q < 10; q++)
    {
      int index = (ola_peak_index - LOW_max_val) + q;
      if( index < 0 )
      {
        index += ola_max;
      } 
      if ( abs( ola_buffer[index] ) > LOW_max_val ) 
      {
        LOW_max_val = abs( ola_buffer[index] );
      }
    }
    
    if( LOW_max_val < (ola_peak/10) )
    {
      detected_thresh += 1;
    }
    
  }
  
}

    
/*MATCHED FILTER*/
void generateTempate(float32_t out_fsk_modulated[][N_LEN], int* in_bit_sequence, int num_bits, float fs, int samples_per_bit, float f0, float f1)
{
  float f_mux;
  float f_chirp = 100*DOWN_SAMP;
  float f_add = 200*DOWN_SAMP;
 
  const int fhss_num = 2;
  float fhss_0[ fhss_num ] = {0};
  float fhss_1[ fhss_num ] = {0};

  int tally_0 = 0;
  int tally_1 = 0;

  for(int m = fhss_num - 1; m >= 0; m--)
  {
    fhss_0[m] = f0 + m*f_add;
    fhss_1[(fhss_num - 1) - m] = f1 - m*f_add;
  }
 
  int modu_len = samples_per_bit*num_bits;
  float32_t pre_modulated[ modu_len ] = {0};
  float32_t reversed_modu[ modu_len ] = {0};
  
  for (int i = 0; i < num_bits; i++)
  {
    int index = (int)(i * samples_per_bit);
    
    if (in_bit_sequence[i] == 0)
    {
      for (int cnt = 0; cnt < samples_per_bit; cnt++)
      {
        f_mux = fhss_0[ tally_0 % fhss_num ] - (f_chirp*cnt/samples_per_bit - f_chirp/2);       
        pre_modulated[index + cnt] = -1*arm_sin_f32( (float32_t)(2 * PI * f_mux * cnt ) / fs );
      }
      tally_0++;
    }
    else
    {
      for (int cnt = 0; cnt < samples_per_bit; cnt++)
      {
        f_mux = fhss_1[ tally_1 % fhss_num ] + (f_chirp*cnt/samples_per_bit - f_chirp/2);
        pre_modulated[index + cnt] = -1*arm_sin_f32( (float32_t)(2 * PI * f_mux * cnt ) / fs );
      }
      tally_1++;
    }

  }
  
  for(int k = 0; k < modu_len; k++)
  {
    reversed_modu[k] = pre_modulated[modu_len - k - 1];
  }

  for (int p = 0; p < B_BLOCKS; p++)
  {
    for (int j = 0; j < L_LEN; j++)
    {
      int indx = p*L_LEN;

      if(indx+j < modu_len)
      {
        out_fsk_modulated[p][j] = reversed_modu[indx + j];
      }
      else
      {
        out_fsk_modulated[p][j] = 0.0;
      }
      
    }
  }
}


/*BLOCK CONVOLVER*/
void UPOLA(float32_t* rx_down)
{
  float32_t y_fft[2*N_LEN] = {0};   
      
  f32tRealFFT(rx_fft[fdata_head], rx_down, N_LEN);

  float32_t sum_msg[N_LEN] = {0};
  float32_t sum_fft[2*N_LEN] = {0};

  
  int index;
  for(int k = 0; k < B_BLOCKS; k++)
  {
    index = fdata_head-k;
    
    if(index < 0)
    {
      index += B_BLOCKS;
    }
    
    arm_cmplx_mult_cmplx_f32(rx_fft[index], mf_fft[k], y_fft, N_LEN);
    
    for(int j = 0; j < 2*N_LEN; j++)
    {
      sum_fft[j] += y_fft[j];
    }
  }
  f32tRealIFFT(sum_msg, sum_fft, N_LEN);



  if(fdata_head >= B_BLOCKS-1)
  {
    fdata_head = 0;
  }
  else
  {
    fdata_head++;
  }

  write_ola_buffer(sum_msg);
}


/*DOWN-SAMPLING*/
void downsample(float32_t* out_arr, float32_t* in_arr, int len)
{
  int cnt = 0;
  for(int k = 0; k < len; k++)
  {
    out_arr[k] = in_arr[cnt];
    cnt += DOWN_SAMP;
  }
}


/*WRITE TO OLA_BUFFER*/
void write_ola_buffer(float32_t* input)
{ 
  ola_head -= ola_overlap;
  
  if(ola_head < 0)
  {
    ola_head += ola_max;
  }
    
  for(int i = 0; i < ola_overlap; i++)
  {
    ola_buffer[ola_head] +=  input[i];

    if ( (abs(ola_buffer[ola_head]) > ola_peak) ) 
    {
      ola_peak_index = ola_head;
      ola_peak = abs(ola_buffer[ola_head]);
    }
    
    ola_head++;
    if(ola_head >= ola_max)
    {
      ola_head = 0;
    }
  }
    
  for(int i = ola_overlap; i < N_LEN; i++)
  {
    ola_buffer[ola_head++] = input[i];
    
    if(ola_head >= ola_max)
    {
      ola_head = 0;
    } 
  }

}


/*FFT f32_t*/
void f32tRealFFT(float32_t* out_fsk_fft, float32_t* in_fsk_msg, int in_len)
{
  arm_rfft_fast_instance_f32 fastfft32;

  arm_rfft_fast_init_f32(&fastfft32, in_len);

  arm_rfft_fast_f32(&fastfft32, in_fsk_msg, out_fsk_fft, 0);
}


/*IFFT f32_t*/
void f32tRealIFFT(float32_t* out_fsk_msg, float32_t* in_fsk_fft, int out_len)
{
  arm_rfft_fast_instance_f32 fastfft32;

  arm_rfft_fast_init_f32(&fastfft32, out_len);

  arm_rfft_fast_f32(&fastfft32, in_fsk_fft, out_fsk_msg, 1);
}
