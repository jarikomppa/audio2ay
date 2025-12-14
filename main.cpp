/* audio 2 ay -- audio files to ay chip commands converter
  version 1.0, December 14th, 2025

  Copyright (C) 2025 Jari Komppa

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
	 claim that you wrote the original software. If you use this software
	 in a product, an acknowledgment in the product documentation would be
	 appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
	 misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jari Komppa
  https://iki.fi/sol

  (i.e, same as zlib license)
*/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <iostream> // for optionparser.h
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fft/soloud_fft.h"
#define DR_WAV_IMPLEMENTATION
#include "dr/dr_wav.h"
#define DR_MP3_IMPLEMENTATION
#include "dr/dr_mp3.h"
#define DR_FLAC_IMPLEMENTATION
#include "dr/dr_flac.h"
#include "stb/stb_vorbis.c"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"
extern "C"
{
#include "ayumi/ayumi.h"
}
#include "optionparser/optionparser.h"

using namespace SoLoud; // for fft

// Options
char* gOutFilename{};
int gFramesPerData{};
int gMinHz{};
int gMaxHz{};
int gAyChannels{};
int gArpLast{};
int gArpRange{};
float gVolBoost{};
int gMonoMinus{};
int gSharpenPasses{};
int gAyRate{};
char* gWavOutFilename{};
int gWavOutMono{};
int gFFTOut{};
char* gWaterfallFilename{};

// Other globals
int gChannels;
int gSampleRate;
int gSamples;
float* gSampleData = NULL;
int gOutSampleCount;
char* gOutSampleData;

unsigned int* gFFTfb{};

// Debug function to render the FFT data and hilight the peaks picked
void render_fft(int block, float* d, int pot, int* hot, int hots)
{
	if (!gFFTfb) gFFTfb = new unsigned int[pot * pot];
	memset(gFFTfb, 0xff, pot * pot * sizeof(unsigned int));
	float maxval = 0;
	for (int i = 0; i < pot; i++)
		if (d[i] > maxval)
			maxval = d[i];
	float mag = 1;
	while (mag < maxval) mag *= 10;

	for (int i = 0; i < pot / 10; i++)
		for (int j = 0; j < pot; j++)
		{
			gFFTfb[(i * 10)*pot + j] = 0xffdddddd;
			gFFTfb[i * 10 + j*pot] = 0xffdddddd;
		}

	for (int i = 0; i < pot / 100; i++)
		for (int j = 0; j < pot; j++)
		{
			gFFTfb[(i * 100) * pot + j] = 0xffaaaaaa;
			gFFTfb[i * 100 + j * pot] = 0xffaaaaaa;
		}

	for (int i = 0; i < pot; i++)
	{
		int ishot = 0;
		for (int j = 0; j < hots; j++)
			if (i == hot[j])
				ishot = 1;
		int ht = (int)(d[i] * pot / mag);
		for (int j = 0; j < ht; j++)
			gFFTfb[i + (pot - j - 1) * pot] = ishot ? 0xff0000ff : 0xff000000;
	}
	char temp[256];
	sprintf(temp, "block%04d.png", block);
	stbi_write_png(temp, pot, pot, 4, gFFTfb, pot * 4);
}

// Convert incoming data to mono.
void convert_to_mono()
{
	if (gChannels == 1) return;

	int multiplier = 1;
	if (gMonoMinus) multiplier = -1;

	for (int i = 0; i < gSamples; i++)
	{
		float t = 0;
		for (int j = 0; j < gChannels; j++)
			t += gSampleData[i * gChannels + j] * multiplier;
		gSampleData[i] = t;
	}
	gChannels = 1;
}

// Yes, it's bubble sort. Not a bottleneck.
void sortby(int *a, float*v, int count)
{
	for (int i = 0; i < count; i++)
	{
		for (int j = i + 1; j < count; j++)
		{
			if (v[a[i]] < v[a[j]])
			{
				int t = a[i];
				a[i] = a[j];
				a[j] = t;
			}
		}
	}
}

// Load wav, ogg, mp3, flac
void load_data(const char* fn)
{
	printf("Loading %s\n", fn);
	FILE* f = fopen(fn, "rb");
	if (!f)
	{
		printf("File not found\n");
		exit(0);
	}
	// Let's assume wav.
	drwav_uint64 totalPCMFrameCount;
	gSampleData = drwav_open_file_and_read_pcm_frames_f32(fn, (unsigned int*)&gChannels, (unsigned int*)&gSampleRate, &totalPCMFrameCount, NULL);

	if (!gSampleData)
	{
		// not a wav, what about mp3?
		drmp3_config conf;
		gSampleData = drmp3_open_file_and_read_pcm_frames_f32(fn, &conf, &totalPCMFrameCount, NULL);
		if (gSampleData)
		{
			gSampleRate = conf.sampleRate;
			gChannels = conf.channels;
		}
	}

	if (!gSampleData)
	{
		// flac, maybe?
		gSampleData = drflac_open_file_and_read_pcm_frames_f32(fn, (unsigned int*)&gChannels, (unsigned int*)&gSampleRate, &totalPCMFrameCount, NULL);
	}

	if (!gSampleData)
	{
		// it must be an ogg then.
		short* output;
		int frames = stb_vorbis_decode_filename(fn, &gChannels, &gSampleRate, &output);
		if (frames > 0)
		{
			unsigned int datalen = frames * gChannels;
			gSampleData = new float[datalen];
			for (unsigned int i = 0; i < datalen; i++)
				gSampleData[i] = output[i] * (1.0f / 0x7fff);
			free(output);
			totalPCMFrameCount = frames;
		}
	}

	if (!gSampleData)
	{
		// okay, keep your secrets.
		printf("Failed to load data\n");
		exit(1);
	}

	gSamples = (unsigned int)totalPCMFrameCount;
}


// Convert 0..1 float volume to 0..15 AY volume
int linear_to_ay_volume(float f)
{
	static const double AY_dac_table[] = {
	  0.0,
	  0.00999465934234,
	  0.0144502937362,
	  0.0210574502174,
	  
	  0.0307011520562,
	  0.0455481803616,
	  0.0644998855573,
	  0.107362478065,
	  
	  0.126588845655,
	  0.20498970016,
	  0.292210269322,
	  0.372838941024,
	  
	  0.492530708782,
	  0.635324635691,
	  0.805584802014,
	  1.0
	};

	int i = 0;
	while (i < 15 && AY_dac_table[i] < f) i++;

	return i;
}

void fputw(int value, FILE* f)
{
	fputc(value & 0xff, f);
	fputc((value >> 8) & 0xff, f);
}

// Actual conversion function
void process()
{
	FILE* f = fopen(gOutFilename, "wb");
	if (!f)
	{
		printf("Unable to open %s\n", gOutFilename);
		exit(1);
	}
	fputc(gAyChannels, f); // ay channels
	fputc(gFramesPerData, f); // output rate

	// Set up AY emulator
	ayumi* ay = new ayumi{};
	ayumi_configure(ay, 0, gAyRate, 44100); // not ym, ay rate (default speccy 1774400), output rate
	if (gWavOutMono)
	{
		ayumi_set_pan(ay, 0, 0.5, 0);
		ayumi_set_pan(ay, 1, 0.5, 0);
		ayumi_set_pan(ay, 2, 0.5, 0);
	}
	else
	{
		ayumi_set_pan(ay, 0, 1,   0);
		ayumi_set_pan(ay, 1, 0.5, 0);
		ayumi_set_pan(ay, 2, 0,   0);
	}
	ayumi_set_mixer(ay, 0, 0, 1, 0); // tone on, noise off, envelope off
	ayumi_set_mixer(ay, 1, 0, 1, 0); // tone on, noise off, envelope off
	ayumi_set_mixer(ay, 2, 0, 1, 0); // tone on, noise off, envelope off

	// Some stats
	int spf = gSampleRate / 50;
	spf *= gFramesPerData;
	int pot = 1;
	while (pot < spf) pot *= 2;
	int blocks = ((gSamples-pot) / spf);
	float* temp = new float[pot*4] {};
	int outspf = 44100 / 50;
	outspf *= gFramesPerData;

	fputw(blocks, f);

	printf("samples per block: %d (at %dhz), power of two window %d\n", spf, gSampleRate, pot);
	printf("samples: %d -> blocks: %d of %d samples\n", gSamples, blocks, outspf);
	printf("in:%3.3fs, out:%3.3fs", gSamples / (float)gSampleRate, blocks * outspf / 44100.0);
	printf(" (little bit will be cropped off the end, it's just how things work)\n");

	// prep output buffers
	gOutSampleCount = (int)(44100.0 * (float)gSamples / (float)gSampleRate);
	gOutSampleData = new char[gOutSampleCount*2];
	memset(gOutSampleData, 127, gOutSampleCount*2);
	float* outtemp = new float[outspf*2];

	float hzpernode = gSampleRate / (float)(pot * 2);
	int minfreq = gMinHz;
	int maxfreq = gMaxHz;
	int minfft = (int)(minfreq / hzpernode);
	int maxfft = (int)(maxfreq / hzpernode);
	if (maxfft > pot) maxfft = pot;
	if (minfft >= maxfft)
	{
		printf("All frequencies filtered out\n");
		exit(0);
	}

	// prep waterfall image
	unsigned int *waterfall = new unsigned int[pot * blocks];
	for (int i = 0; i < pot * blocks; i++)
		waterfall[i] = 0xff777777;

	printf("Using band %d-%dhz\n", (int)(minfft * hzpernode), (int)(maxfft * hzpernode));

	int* best = new int[pot]{};
	int lowhz = 1000000;
	int highhz = 0;
	float global_peak = 0;
	// need two passes to find out how loud we are for output (second pass does actual processing)
	for (int pass = 0; pass < 2; pass++)
	for (int i = 0; i < blocks; i++)
	{
		// perform stft
		memset(temp, 0, sizeof(float) * pot * 4);
		for (int j = 0; j < pot; j++)
		{
			temp[j*2] = gSampleData[spf * i + j];
		}
		FFT::fft(temp, pot*2);
		for (int j = 0; j < pot; j++)
		{
			float real = temp[j * 2];
			float imag = temp[j * 2 + 1];
			temp[j] = (float)sqrt(real * real + imag * imag);
		}

		// find the global peak on first pass
		if (pass == 0)
		{
			for (int j = minfft; j < maxfft-minfft; j++)
				if (global_peak < temp[j])
					global_peak = temp[j];
		}

		// do the actual work on second pass
		if (pass == 1)
		{
			memcpy(temp + pot, temp, pot * sizeof(float));
			
			// sharpen filter
			for (int k = 0; k < gSharpenPasses; k++)
			{
				memcpy(temp + pot * 2, temp + pot, pot * sizeof(float));
				temp[pot + 0] = 0;
				temp[pot + 1] = 0;
				temp[pot + pot - 1] = 0;
				temp[pot + pot - 2] = 0;
				for (int j = 2; j < pot - 2; j++)
				{
					temp[pot + j] = temp[pot * 2 + j] * 4 - 
						            temp[pot * 2 + j - 1] - 
						            temp[pot * 2 + j - 2] - 
						            temp[pot * 2 + j + 1] - 
						            temp[pot * 2 + j + 2];
					if (temp[pot + j] < 0) temp[pot + j] = 0;
				}
			}
			// temp+pot = sharpen filtered
				
			for (int j = 0; j < outspf*2; j++)
				outtemp[j] = 0;
			for (int j = 0; j < pot; j++)
				best[j] = j + minfft;
			sortby(best, temp + pot, maxfft - minfft); // use sharpened values for finding peaks

			// default volumes to zero in case all channels are not in use
			ayumi_set_volume(ay, 0, 0);
			ayumi_set_volume(ay, 1, 0);
			ayumi_set_volume(ay, 2, 0);

			// and now, for three best peaks (more or less), use AY for reproduction
			int nextpeak = 0;
			for (int j = 0; j < gAyChannels; j++)
			{
				int peak = best[nextpeak];
				int skips = 1;				

				if (gArpLast && j == gAyChannels-1) skips = i % gArpRange;
				for (int k = 0; k < skips; k++)
				{
					while (nextpeak > 0 && abs(best[nextpeak] - best[nextpeak - 1]) == 1) nextpeak++;
					nextpeak++;
				}
				
				if (gAyChannels == 1 && gArpLast) // special case: first is last and we want to arp
				{
					peak = nextpeak - 1; // not exactly the same but close enough
					if (peak < 0) peak = 0;
					peak = best[peak];
				}

				float hz = 0;
				float hzd = 0;
				float volume = 0;
				// use weighted average of nearby values to find the "exact" frequency; this merges peaks
				for (int k = 0; k < 5; k++)
				{
					if (peak >= 2 && peak < pot - 2)
					{
						hz += (peak - 2 + k) * temp[(peak - 2 + k)];
						hzd += temp[(peak - 2 + k)];
					}
				}
				if (hzd > 0)
				{
					hz = hz * hzpernode / hzd;
				}

				if (hz < minfreq)
				{
					volume = 0;
				}
				else
				{
					volume = (hzd / 2) / global_peak;
					if (hz > highhz)
						highhz = (int)hz;
					if (hz < lowhz)
						lowhz = (int)hz;
				}
				
				volume *= gVolBoost;

				int intvol = linear_to_ay_volume(volume);

				int intensity = (int)(0xff * volume);
				waterfall[blocks * (pot-best[j]) + i] = 0xff000000 + 0x01 * intensity;
				
				int tone = (gAyRate / 16) / (int)hz;
				
				unsigned short outvalue = (tone & 0xfff) | (intvol << 12);
				fputw(outvalue, f);


				ayumi_set_volume(ay, j, (outvalue >> 12) & 0xf);
				ayumi_set_tone(ay, j, outvalue & 0xfff);
			}

			// render audio
			for (int k = 0; k < outspf; k++)
			{
				ayumi_process(ay);
				outtemp[k * 2 + 0] = (float)ay->left / 3;
				outtemp[k * 2 + 1] = (float)ay->right / 3;
			}
			// and from float audio to 8 bit
			for (int j = 0; j < outspf*2; j++)
			{
				gOutSampleData[i * outspf*2 + j] = 127 + (signed char)(outtemp[j] * 100);
			}
		}
		if (gFFTOut)
			render_fft(i, temp, pot, best, 3);
	}
	printf("Ended up using: %d-%dhz\n", lowhz, highhz);
	if (gWaterfallFilename != 0)
		stbi_write_png(gWaterfallFilename, blocks, pot, 4, waterfall, blocks * 4);
	fclose(f);
	delete[] waterfall;
	delete[] outtemp;
	delete ay;
	delete[] temp;
	delete[] best;
}

// Save resulting simulated audio
void save_wave()
{
	if (gWavOutFilename == 0)
		return;
	drwav_data_format format;
	format.container = drwav_container_riff;     // <-- drwav_container_riff = normal WAV files, drwav_container_w64 = Sony Wave64.
	format.format = DR_WAVE_FORMAT_PCM;          // <-- Any of the DR_WAVE_FORMAT_* codes.
	format.channels = 2;
	format.sampleRate = 44100;
	format.bitsPerSample = 8;
	drwav wav;	
	printf("saving %s\n", gWavOutFilename);
	if (!drwav_init_file_write(&wav, gWavOutFilename, &format, NULL))
	{
		printf("Unable to open output file\n");
		exit(0);
	}
	drwav_write_pcm_frames(&wav, gOutSampleCount, gOutSampleData);
	drwav_uninit(&wav);
}

enum optionIndex { UNKNOWN, HELP, OUTFILENAME, FRAMESPERDATA, MINHZ, MAXHZ, AYCHANNELS, ARPLAST, ARPRANGE, VOLBOOST, MONOMINUS, SHARPENPASSES, AYRATE, WAVOUT, WAVOUTMONO, FFTOUT, WATERFALLOUT, SHOWOPT};
const option::Descriptor usage[] =
{
	{ UNKNOWN,		0, "", "",	option::Arg::None,				 "USAGE: audio2ay inputfilename [options]\n\nOptions:"},
	{ HELP,			0, "h", "help", option::Arg::None,			 " -h --help\t Print usage and exit"},
	{ OUTFILENAME,  0, "o", "output", option::Arg::Optional,     " -o --output=filename\t output filename, default \"aydata.dat\""},
	{ FRAMESPERDATA,0, "f", "frames", option::Arg::Optional,     " -f --frames=num\t Frames per data, default 1"},
	{ MINHZ,        0, "a", "minhz", option::Arg::Optional,      " -a --minhz=num\t Minimum frequency to consider, default 300"},
	{ MAXHZ,        0, "b", "maxhz", option::Arg::Optional,      " -b --maxhz=num\t Maximum frequency to consider, default 2000"},
	{ AYCHANNELS,   0, "c", "channels", option::Arg::Optional,   " -c --channels=num\t AY channels to use, default 3"},
	{ ARPLAST,      0, "p", "arplast", option::Arg::Optional,    " -p --arplast=num\t Arpeggiate last channel, default 1 if more than 1 channel"},
	{ ARPRANGE,     0, "g", "arprange", option::Arg::Optional,   " -g --arprange=num\t Peaks to cycle, default 3"},
	{ VOLBOOST,     0, "v", "volboost", option::Arg::Optional,   " -v --volboost=num\t Volume multiplier, default 1.5"},
	{ MONOMINUS,    0, "m", "monominus", option::Arg::Optional,  " -m --monominus=num\t Use - instead of + when combining channels, default 0"},
	{ SHARPENPASSES,0, "s", "sharpen", option::Arg::Optional,    " -s --sharpen=num\t Sharpen filter passes on FFT data, default 4"},
	{ AYRATE,       0, "y", "ayrate", option::Arg::Optional,     " -y --ayrate=num\t AY chip rate, default 1774400 (zx spectrum)"},
	{ WAVOUT,       0, "w", "wavout", option::Arg::Optional,     " -w --wavout=filename\t Output simulated wav file, default none"},
	{ WAVOUTMONO,   0, "z", "wavoutmono", option::Arg::Optional, " -z --wavoutmono=num\t Use mono panning on wav output, default=0"},
	{ FFTOUT,       0, "t", "fftout", option::Arg::Optional,     " -t --fftout=num\t Output tons of fft images, default=0"},
	{ WATERFALLOUT, 0, "l", "waterfall", option::Arg::Optional,  " -l --waterfall=filename\t Output waterfall image on chosen notes, default=0"},
	{ SHOWOPT,      0, "q", "showopt", option::Arg::Optional,    " -q --showopt=num\t Show used options, default=0"},
	{ UNKNOWN,      0, "", "", option::Arg::None,				 "Example:\n  audio2ay test.wav -o=test.dat -w=test.wav\nSupported input file types: wav, flac, mp3, ogg\n"},
	{ 0,0,0,0,0,0 }
};


int main(int parc, char** pars)
{
	printf("Audio to AY encoder by Jari Komppa 2025 http://iki.fi/sol\n\n");

	option::Stats stats(usage, parc - 1, pars + 1);
	assert(stats.buffer_max < 32 && stats.options_max < 32);
	option::Option options[32], buffer[32];
	option::Parser parse(true, usage, parc - 1, pars + 1, options, buffer);

	if (options[UNKNOWN])
	{
		for (option::Option* opt = options[UNKNOWN]; opt; opt = opt->next())
			printf("Unknown option: %s\n", opt->name);
		printf("Run without parameters for help.\n");
		exit(0);
	}

	if (parse.error() || parc < 2 || options[HELP] || parse.nonOptionsCount() != 1)
	{
		option::printUsage(std::cout, usage);
		return 0;
	}

#define INTPARAM(ENUM,VAR,DEFAULT,MIN,MAX)	\
	if (options[ENUM] && options[ENUM].arg)	\
	{										\
		VAR = atoi(options[ENUM].arg);		\
		if (VAR > MAX || VAR < MIN)			\
		{									\
			printf("Invalid value for %s: %d (%d - %d)\n", #ENUM,VAR, MIN, MAX); \
			exit(1);						\
		}									\
	}										\
	else									\
	{										\
		VAR = DEFAULT;						\
	}
#define FLOATPARAM(ENUM,VAR,DEFAULT,MIN,MAX)		\
	if (options[ENUM] && options[ENUM].arg)	\
	{										\
		VAR = (float)atof(options[ENUM].arg);		\
		if (VAR > MAX || VAR < MIN)			\
		{									\
			printf("Invalid value for %s: %3.3f (%3.3f-%3.3f)\n", #ENUM,VAR,(float)MIN,(float)MAX); \
			exit(1);						\
		}									\
	}										\
	else									\
	{										\
		VAR = DEFAULT;						\
	}
#define STRINGPARAM(ENUM,VAR,DEFAULT)		\
	if (options[ENUM] && options[ENUM].arg)	\
	{										\
		VAR = _strdup(options[ENUM].arg);	\
	}										\
	else									\
	{										\
		VAR = (char*)DEFAULT;				\
	}

	STRINGPARAM(OUTFILENAME, gOutFilename, "aydata.dat")
	INTPARAM(FRAMESPERDATA, gFramesPerData, 1, 1, 100)
	INTPARAM(MINHZ, gMinHz, 300, 0, 44100)
	INTPARAM(MAXHZ, gMaxHz, 2000, 0, 44100)
	INTPARAM(AYCHANNELS, gAyChannels, 3, 1, 3)
	INTPARAM(ARPLAST, gArpLast, (gAyChannels == 1 ? 0 : 1), 0, 1)
	INTPARAM(ARPRANGE, gArpRange, 3, 1, 100)
	FLOATPARAM(VOLBOOST, gVolBoost, 1.5, 0, 100)
	INTPARAM(MONOMINUS, gMonoMinus, 0, 0, 1)
	INTPARAM(SHARPENPASSES, gSharpenPasses, 4, 0, 100)
	INTPARAM(AYRATE, gAyRate, 1774400, 100000, 10000000)
	STRINGPARAM(WAVOUT, gWavOutFilename, 0)
	INTPARAM(WAVOUTMONO, gWavOutMono, 0, 0, 1)
	INTPARAM(FFTOUT, gFFTOut, 0, 0, 1)
	STRINGPARAM(WATERFALLOUT, gWaterfallFilename, 0)
	int showopt = 0;
	INTPARAM(SHOWOPT,showopt,0,0,1)

	if (showopt)
	{
		printf("Options in use:\n"
			"Output filename: %s\n"
			"Frames per data: %d\n"
			"MinHz          : %d\n"
			"MaxHz          : %d\n"
			"AY Channels    : %d\n"
			"Arp last       : %d\n"
			"Arp range      : %d\n"
			"Volume boost   : %3.3f\n"
			"Mono-minus     : %d\n"
			"Sharp passes   : %d\n"
			"AY rate        : %d\n"
			"Wave out       : %s\n"
			"Wave out mono  : %d\n"
			"FFT output     : %d\n"
			"Waterfall      : %s\n\n",
			gOutFilename,
			gFramesPerData,
			gMinHz,
			gMaxHz,
			gAyChannels,
			gArpLast,
			gArpRange,
			gVolBoost,
			gMonoMinus,
			gSharpenPasses,
			gAyRate,
			gWavOutFilename ? gWaterfallFilename : "< no wave output >",
			gWavOutMono,
			gFFTOut,
			gWaterfallFilename ? gWaterfallFilename : "< no waterfall output >");
	}

	load_data(parse.nonOption(0));
	convert_to_mono();
	process();
	save_wave();
	printf("all done, %s written\n", gOutFilename);
	return 0;
}