# Speech coursework 1

This is my code for the first coursework in EEEM030 (Speech and Audio Processing) at the University of Surrey. It is a vowel synthesiser. It takes the vowel in "who'd" from a male and a female recording, works out the shape of the vocal tract from 100 ms of each with linear predictive coding (LPC), and then builds a one-second synthetic vowel by passing a train of pulses at the same pitch through that filter.

I also measure the first three formants and the pitch of each vowel, and look at how the result changes with the LPC order and the length of the segment. The report is in the docs folder and explains all of this.

## Running it

It is written in Matlab (R2016b or later, and plain Matlab is enough). Open `cw1_lpc_vowels.m` and press Run. It takes about a minute and puts everything in the `output` folder, replacing what is already there. If the Signal Processing Toolbox is installed, it also compares my LPC code with the built-in `lpc` and `arcov`, but it does not need them.

You can also run it online without installing anything. Open this link and sign in with your university account:
https://matlab.mathworks.com/open/github/v1?repo=TheFinix13/speech-coursework-1

When it finishes, the last lines say "6 of 6 self-tests passed" and "Done". For the male vowel I get formants of 313, 1230 and 2264 Hz and a pitch of 124.5 Hz. For the female vowel I get 382, 1842 and 2378 Hz and 200.0 Hz.

## What is where

- `cw1_lpc_vowels.m` is the whole analysis, with the functions at the bottom of the file.
- `data/speech` holds the recordings I used, from the sample set for the assignment. `hod_f.wav` was my first choice for the female vowel, and I kept it as an extra.
- `output/figures` holds the figures in the report.
- `output/tables` holds the tables as csv files, all the numbers in `results.json`, and `summary.txt`, which is what Matlab printed when I ran it.
- `output/audio` holds the synthesised vowels. `whod_m_synth_final.wav` and `whod_f_synth_final.wav` are the main two. The other files change one thing: `order04` to `order40` use a different LPC order, `len020ms` to `len200ms` use a shorter or longer segment, `preemph` adds pre-emphasis, and `10kHz_order12` resamples to 10 kHz first. `original_segment` is the recorded piece, and `impulse_train` is the pulses before filtering.
- `docs` holds the report as a PDF.

Fiyinfoluwa Akano, MSc Artificial Intelligence, University of Surrey
