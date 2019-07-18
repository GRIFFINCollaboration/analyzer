/******************************************************************************
 *       Simple DESCANT analyser
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h> /* rand() */
#include <time.h>
#include <unistd.h> /* usleep */
#include <math.h>
#include <string.h>  /* memset */
#include <stdint.h>
#include <stdbool.h>

#include "midas.h"
#include "histogram.h"
#include "web_server.h"

#define TRUE 1
#define FALSE 0

int descant(EVENT_HEADER *, void *);
int descant_init();
int descant_bor(INT run_number);
int descant_eor(INT run_number);
extern void read_odb_gains();
extern int read_odb_histinfo();
extern int hist_init();
extern float spread(int val);

ANA_MODULE descant_module = {
	"Descant",              /* module name           */
	"UrMom",                /* author                */
	descant,                /* event routine         */
	descant_bor,            /* BOR routine           */
	descant_eor,            /* EOR routine           */
	descant_init,           /* init routine          */
	NULL,                   /* exit routine          */
	NULL,                   /* parameter structure   */
	0,                      /* structure size        */
	NULL,                   /* initial parameters    */
};

#define NUM_ODB_CHAN     459 // size of msc table in odb
#define MAX_SAMPLE_LEN  4096
#define ENERGY_BINS    65536 /* 65536 131072 262144 */
#define NUM_CLOVER        16
#define MAX_CHAN        1024 
#define STRING_LEN       256
#define MIDAS_STRLEN      32
#define MAX_ADDRESS  0x10000
#define E_SPEC_LENGTH   8192
#define T_SPEC_LENGTH   8192
#define WV_SPEC_LENGTH  4096

typedef struct descant_fragment_struct {
	int   address;
	int   chan;
	long  timestamp;
	int   cfd;
	int   energy;
	int   overrange;
	int   cc_short;
	int   tot_trig;
	int   lost_trig;
	short waveform_length;  
	int   flags;
} Desc_event;

static Desc_event  desc_event;
static short       waveform1[MAX_SAMPLE_LEN];
static short       waveform2[MAX_SAMPLE_LEN];
static short       digitalWaveform1[MAX_SAMPLE_LEN];
static short       digitalWaveform2[MAX_SAMPLE_LEN];
static int        rate_data[MAX_CHAN];
static int     num_chanhist;

extern HNDLE hDB; // Odb Handle

float   gains[NUM_ODB_CHAN];
float offsets[NUM_ODB_CHAN];

int decode_descant_event( unsigned int *evntbuf, int evntbuflen);
int process_decoded_descant(Desc_event *ptr);
int unpack_descant_bank(unsigned *buf, int len);

/////////////////////////////////////////////////////////////////////////////

int descant_init()
{  int i;
	read_odb_gains();  // Print the loaded gains and offsets
	fprintf(stdout,"\nRead Gain/Offset values from ODB\nIndex\tGain\tOffset\n");
	for(i=0; i<NUM_ODB_CHAN; i++){
		fprintf(stdout,"%d\t%f\t%f\n",i,gains[i],offsets[i]);
	}
	fprintf(stdout,"\n\n");

	if( (num_chanhist = read_odb_histinfo()) <= 0 ){
		;
	}
	hist_init();
	return SUCCESS;
}
int descant_eor(INT run_number){ return SUCCESS; }

extern TH1I **hit_hist;
extern TH1I **sum_hist;
extern TH1I **ph_hist;
extern TH1I **e_hist;
extern TH1I **cfd_hist;
extern TH1I **wave_hist;

static time_t bor_time=-1, start_time=-1;

int descant_bor(INT run_number)
{
	bor_time = time(NULL); start_time = -1;
	read_odb_gains();
	Zero_Histograms();
	memset(rate_data, 0, sizeof(rate_data));
	printf("Success!\n");
	return SUCCESS;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////            per event functions           //////////////////
/////////////////////////////////////////////////////////////////////////////

int descant(EVENT_HEADER *pheader, void *pevent)
{
	unsigned int *data;
	int words;

	if(        (words = bk_locate(pevent, "GRF0", (DWORD *)&data)) > 0 ){
	} else if( (words = bk_locate(pevent, "GRF3", (DWORD *)&data)) > 0 ){
	} else if( (words = bk_locate(pevent, "GRF4", (DWORD *)&data)) > 0 ){
	} else {
		return(-1);
	}
	decode_descant_event(data, words);
	return SUCCESS;
}

extern int GetIDfromAddress(int addr);

int decode_descant_event(uint32_t* data, int size)
{
	int w = 0;
	int board;
	uint8_t channel;
	int ev;
	int s;

	// Clear desc_event at start of buffer and check the first word is a good header
	memset(&desc_event, 0, sizeof(Desc_event) );

	for(board = 0; w < size; ++board) {
		// read board aggregate header
		if(data[w]>>28 != 0xa) {
			return 0;
		}
		int32_t numWordsBoard = data[w++]&0xfffffff; // this is the number of 32-bit words from this board
		if(w - 1 + numWordsBoard > size) {
			return 0;
		}
		uint8_t boardId = data[w]>>27; // GEO address of board (can be set via register 0xef08 for VME)
		//uint16_t pattern = (data[w]>>8) & 0x7fff; // value read from LVDS I/O (VME only)
		uint8_t channelMask = data[w++]&0xff; // which channels are in this board aggregate
		++w;//uint32_t boardCounter = data[w++]&0x7fffff; // ??? "counts the board aggregate"
		++w;//uint32_t boardTime = data[w++]; // time of creation of aggregate (does not correspond to a physical quantity)
		// Loop over all even channels (odd channels are grouped with the corresponding even channel)
		for(channel = 0; channel < 16; channel += 2) {
			// skip channels not in the channel mask
			if(((channelMask>>(channel/2)) & 0x1) == 0x0) {
				continue;
			}
			// read channel aggregate header
			if(data[w]>>31 != 0x1) {
				// bad header, skip rest of data
				return 0;
			}
			int32_t numWords = data[w++]&0x3fffff;//per channel
			if(w >= size) {
				return 0;
			}
			if(((data[w]>>29) & 0x3) != 0x3) {
				return 0;
			}

			// we don't care about dual trace, everything is written to one waveform anyways
			bool dualTrace        = ((data[w]>>31) == 0x1);
			bool extras           = (((data[w]>>28) & 0x1) == 0x1);
			bool waveformPresent  = (((data[w]>>27) & 0x1) == 0x1);
			uint8_t extraFormat   = ((data[w]>>24) & 0x7);
			//for now we ignore the information which traces are stored:
			//bits 22,23: if(dualTrace) 00 = "Input and baseline", 01 = "CFD and Baseline", 10 = "Input and CFD"
			//            else          00 = "Input", 01 = "CFD"
			//bits 19,20,21: 000 = "Long gate",  001 = "over thres.", 010 = "shaped TRG", 011 = "TRG Val. Accept. Win.", 100 = "Pile Up", 101 = "Coincidence", 110 = reserved, 111 = "Trigger"
			//bits 16,17,18: 000 = "Short gate", 001 = "over thres.", 010 = "TRG valid.", 011 = "TRG HoldOff",           100 = "Pile Up", 101 = "Coincidence", 110 = reserved, 111 = "Trigger"
			int numSampleWords = 4*(data[w++]&0xffff);// this is actually the number of samples divided by eight, 2 sample per word => 4*
			if(w >= size) {
				return 0;
			}
			int eventSize = numSampleWords+2; // +2 = trigger time words and charge word
			if(extras) ++eventSize;
			if(numWords%eventSize != 2 && !(eventSize == 2 && numWords%eventSize == 0)) { // 2 header words plus n*eventSize should make up one channel aggregate
				return 0;
			}
			if(dualTrace) desc_event.waveform_length =   numSampleWords;
			else          desc_event.waveform_length = 2*numSampleWords;

			// read channel data
			for(ev = 0; ev < (numWords-2)/eventSize; ++ev) { // -2 = 2 header words for channel aggregate
				//event->SetDaqTimeStamp(boardTime);
				desc_event.address = (0x8000 + (boardId * 0x100) + channel + (data[w]>>31)); // highest bit indicates odd channel
				desc_event.chan = GetIDfromAddress(desc_event.address);
				// these timestamps are in 2ns units
				desc_event.timestamp = (data[w] & 0x7fffffff);
				++w;
				if(waveformPresent) {
					if(w + numSampleWords >= size) { // need to read at least the sample words plus the charge/extra word
						return 0;
					}
					for(s = 0; s < numSampleWords && w < size; ++s, ++w) {
						if(dualTrace) {
							// 14 low bits are first waveform, the 14 high bits the second waveform
							waveform1[s] = data[w]       & 0x3fff;
							waveform2[s] = (data[w]>>16) & 0x3fff;
							// bits 15/16 are first sample of digital waveforms 1/2, bit 30/31 the second
							digitalWaveform1[2*s]   = (data[w]>>14) & 0x1;
							digitalWaveform2[2*s]   = (data[w]>>15) & 0x1;
							digitalWaveform1[2*s+1] = (data[w]>>30) & 0x1;
							digitalWaveform2[2*s+1] = (data[w]>>31) & 0x1;
						} else {
							// 14 low bits are first sample, the 14 high bits the second sample
							waveform1[2*s]   = data[w]       & 0x3fff;
							waveform1[2*s+1] = (data[w]>>16) & 0x3fff;
							// bits 15/16 are first sample of digital waveforms 1/2, bit 30/31 the second
							digitalWaveform1[2*s]   = (data[w]>>14) & 0x1;
							digitalWaveform2[2*s]   = (data[w]>>15) & 0x1;
							digitalWaveform1[2*s+1] = (data[w]>>30) & 0x1;
							digitalWaveform2[2*s+1] = (data[w]>>31) & 0x1;
						}
					}
				} else {
					if(w >= size) { // need to read at least the sample words plus the charge/extra word
						return 0;
					}
				}
				if(extras) {
					switch(extraFormat) {
						case 0: // [31:16] extended time stamp, [15:0] baseline*4
							desc_event.timestamp |= ((uint64_t)data[w]&0xffff0000)<<15;
							break;
						case 1: // [31:16] extended time stamp, 15 trigger lost, 14 over range, 13 1024 triggers, 12 n lost triggers
							desc_event.timestamp |= ((uint64_t)data[w]&0xffff0000)<<15;
							desc_event.flags = (data[w]>>12)&0xf;
							break;
						case 2: // [31:16] extended time stamp,  15 trigger lost, 14 over range, 13 1024 triggers, 12 n lost triggers, [9:0] fine time stamp
							desc_event.timestamp |= ((uint64_t)data[w]&0xffff0000)<<15;
							desc_event.cfd = data[w]&0x3ff;
							desc_event.flags = (data[w]>>12)&0xf;
							break;
						case 4: // [31:16] lost trigger counter, [15:0] total trigger counter
							desc_event.lost_trig = data[w]&0xffff;
							desc_event.tot_trig = data[w]>>16;
							break;
						case 5: // [31:16] CFD sample after zero cross., [15:0] CFD sample before zero cross.
							w++;
							break;
						case 7: // fixed value of 0x12345678
							if(data[w] != 0x12345678) {
								fprintf(stderr,"Failed to get debug data word 0x12345678, got 0x%08x\n", data[w]);
								break;
							}
							break;
						default:
							break;
					}
					++w;
				}
				desc_event.cc_short = data[w]&0x7fff;
				desc_event.overrange = (data[w]>>15) & 0x1;//this is actually the over-range bit!
				desc_event.energy = data[w++]>>16;
				if(!(process_decoded_descant(&desc_event))) return -1;
				memset(&desc_event, 0, sizeof(Desc_event) );
			} // while(w < size)
		} // for(uint8_t channel = 0; channel < 16; channel += 2)
	} // for(int board = 0; w < size; ++board)

	return(SUCCESS); 
}

#define RESET_SECONDS     10
static time_t last_reset;
#define TEST_HITPAT(var,bit) (var[(int)bit/32] & (1<<(bit%32)))
int process_decoded_descant(Desc_event *ptr)
{
	time_t current_time = time(NULL);
	static int last_sample, event;
	int i, chan, ecal, len;
	float energy;
	short sample;

	if( bor_time == -1 ){ bor_time = current_time; }
	if( start_time == -1 ){
		start_time=current_time; event=0;
	} else { ++event; }

	// this should be done in a scalar routine
	if( current_time - last_reset > RESET_SECONDS ){
		memset(rate_data, 0, sizeof(rate_data));
		last_reset = current_time;
	}

	if( (chan = ptr->chan) == -1 ){ return(0); } // msg printed above for these
	if( chan >= num_chanhist ){
		fprintf(stderr,"process_event: ignored event in chan:%d [0x%04x]\n",
				chan, ptr->address );
		return(0);
	}
	energy = ptr->energy;
	ecal   = spread(energy) * gains[chan] + offsets[chan];

	++rate_data[chan];
	ph_hist[chan] -> Fill(ph_hist[chan],  (int)energy,     1);
	hit_hist[3]   -> Fill(hit_hist[3],    chan,            1);
	e_hist[chan]  -> Fill(e_hist[chan],   (int)ecal,       1);
	hit_hist[0]   -> Fill(hit_hist[0],    chan,            1);
	if(chan<64){sum_hist[0]   -> Fill(sum_hist[0],    (int)ecal,       1);} //Quick hack for only LO gain channels
	if(chan>63 && chan<128){sum_hist[1]   -> Fill(sum_hist[1],    (int)ecal,       1);} //Quick hack for only HI gain channels
	if(chan>99 && chan<106){sum_hist[3]   -> Fill(sum_hist[3],    (int)ecal,       1);} //Quick hack for only PACES channels
	if(chan>91 && chan<100){sum_hist[4]   -> Fill(sum_hist[4],    (int)ecal,       1);} //Quick hack for only LaBr3 channels
	cfd_hist[chan]-> Fill(cfd_hist[chan], ptr->cfd/256,    1);
	hit_hist[1]   -> Fill(hit_hist[1],    chan,            1);


	if( (len = ptr->waveform_length) != 0 ){
		if( len > MAX_SAMPLE_LEN ){ len = MAX_SAMPLE_LEN; }

		if( last_sample == 0 ){ last_sample = len; }
		else if( last_sample != len ){
			fprintf(stderr,"event %4d - samples changed %d to %d\n",
					event, last_sample, len);
		}

		hit_hist[2]   -> Fill(hit_hist[2],    chan,            1);
		wave_hist[chan]->Reset(wave_hist[chan]);
		wave_hist[chan]->SetValidLen(wave_hist[chan], len);
		for(i=0; i<len; i++){
			sample = waveform1[i];
			wave_hist[chan]->SetBinContent(wave_hist[chan], i, sample );
		}
	}

	return(0);
}
