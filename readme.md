# Audio to AY

This is a tool that converts any audio stream into low-rate AY-3-8910 data stream.

It does not play the sound as samples, but instead performs short-time fourier transform on short clips of audio, finds the strongest signals and tries to reproduce them with the sound chip, usually badly.

Melodic sources work the best. Speech, for instance, doesn't really work, but maybe the result is what you're after.

## License

The tool itself is licensed under zlib/libpng license:

-- 8< -- 8< -- 8< --

audio 2 ay -- audio files to ay chip commands converter
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

-- 8< -- 8< -- 8< --

The tool uses various other libraries, including dr and stb ones, with
their own licenses. Look into the source files for details. Everything
is pretty liberally licensed, though.

## File Format

1 byte      - number of AY channels used (1-3)
1 byte      - frames per data block (usually 1)
2 bytes     - blocks in the file
2 bytes x N - AY register data. 4 bits volume, 4 bits coarse, 8 bits fine

## Tools for Better Results

To hear the results without going through actual hardware, you can tell
the tool to generate wave files directly.

If the result sounds like random noise, ou can try to hone the results 
by limiting the frequency range. You can also try changing the volume
boost and sharpen filters.

You can generate a waterfall image to see what notes were picked, and
you can also export FFT images per block. Do note that this will generate
very, very many images.

To generate smaller data, you can increase the frames per block count,
and/or reduce the AY channel count.

## Commandline Options

audio2ay inputfilename (and optionally options)

The input file can be a wav, flac, mp3 or ogg file.

Example:

    audio2ay test.wav --output=test.dat --wavout==test.wav


--help               Print usage and exit

Prints out the usage, which is more or less this list, without the
additional explanations.

--output=filename    output filename, default "aydata.dat"

Filename for the AY binary. Always written, if not specified,
uses the default filename.

--frames=num         Frames per data, default 1

Update the AY registers every N frames. Bigger value means smaller
output file, but worse representation.

--minhz=num          Minimum frequency to consider, default 300
--maxhz=num          Maximum frequency to consider, default 2000

These two can be used to limit the frequency ranges to consider. Most sounds
have very loud signals in the bottom range, so it makes sense to remove them.
The top end is generally just noise.

--channels=num       AY channels to use, default 3

By default, use all three channels. 2 channels also work surprisingly well 
(when it works).

--arplast=num        Arpeggiate last channel, default 1 if more than 1 channel

Arpeggiate the last channel. Instead of using the best peak, skip peaks in a
round-robin way (like +0, +1, +2, +0, +1, +2...)

If only one AY channel is used, arpeggiation is disabled by default, otherwise
it's enabled by default.

--arprange=num       Peaks to cycle, default 3

Arpeggiation range, i.e, how how many peaks to skip at maximum.

--volboost=num       Volume multiplier, default 1.5

How much to boost the volume before converting to AY volume. If your result
is way too loud, scale this down.

--monominus=num      Use - instead of + when combining channels, default 0

Do left minus right when converting input data to mono instead of left plus right.

--sharpen=num        Sharpen filter passes on FFT data, default 4

How many passes of sharpen filter should be applied to the FFT data before selecting
the peaks. Removes noise.

--ayrate=num         AY chip rate, default 1774400 (zx spectrum)

AY chip rate to use in calculations. If you're using this to play back on some
other AY chip containing computer, like the Atari ST, you need to give a different
frequency value.

--wavout=filename    Output simulated wav file, default none

Filename for the simulated audio preview. By default, not generated, but useful when
finding nice options.

--wavoutmono=num     Use mono panning on wav output, default=0

Generate mono panning instead of left/center/right panning for the channels in the preview.

--fftout=num         Output tons of fft images, default=0

Enable generation of crapton of .png files, one for each block, that has a drawing of
the FFT data the code is looking at. You probably don't want this.

--waterfall=filename Output waterfall image on chosen notes, default=0

Enable generation of a waterfall image that shows which frequencies were picked. May
be helpful in debugging. Or may not.

--showopt=num        Show used options, default=0

If you're unsure whether the tool is reading your options right, this can be used to
show what it went with in the end.

## ZX Spectrum Example

The example found in the z80 directory compiles into ZX Spectum .tap file.

The example includes the binary file, writes to the AY registers and then
waits for N frames to pass, and repeats. When the data reaches end, it
starts over from the beginning.

Typically you would do the writes in a timer interrupt.

## EOF
