# HMM POS tagger with suffix analysis
### by Santichai Pornavalai
1.4.19

This is a simple program that only takes two commandline arguments
*hmm <PATH TO TRAIN FILE>  <PATH TO TEST FILE>* 

The training file should be in the form of 
    <TOKEN>   <TAG>     this should be separated by exactly 3 spaces... so sorry. 

make tests trains and tags the testing set returns the accuracy. If scikit-learn is installed, a report showing
detailed precision recall values for individual tags would be printed out


Training and decoding is done by maximum likelihood estimation and viterbi algorithm respectively. 
Tag bigrams are smoothed by interpolating with the unigram probabilites (default set to 0.95 for bigrams).
Unknown words are handled using suffix analysis similar to that of TnT. I however tried to work around this by using 
Only one Trie to and one Array for all the suffixes from length n to 0. To combat the frequency imbalance I simply use a 
pseudo count of a number less than 1 instead. I'm not sure if this works but it sums up to 1. The main problem I found
was that, since the program learns suffix probabilites in the same corpus pass as the other probabilites, I can't filter out
which words to perform suffix analysis on before hand. This leads to suffixes of frequent words having very high probabilities. 
Another caveat is that I am basically using Viterbi algorithm with two probabilites, namely emmision and suffix. I found it 
difficult to balance out the two probabilities when combining them with simple interpolation. suffix probabilities are just 
way higher than emmision. 

Indices for words and tags are stored in a Trie modified to handle suffixes as well. Integer to word mappings are stored in 
a character array separated by NULL to minimize space. 

Performance wise it takes around 1-2 seconds to train on the tiger corpus. Tagging isn't too slow with around 6000 - 1000 words tagged per second. This depends on the number of unknown words.
Space complexity isn't too bad either but might be a bit larger than usual because I'm only using double floats to be safe.

I didn't perform any rigorous performance testing but splitting the tiger corpus to a testing and training corpus (1/5 split) yields 89% accuracy and similar F1 score.
This is not on par with the original TnT tagger. The reason maybe: my (mis)interpretation of suffix handling, only bigram tags (TNT uses trigrams), different training corpus etc.


There is a bug, in which the tag probabilities of the punctuations don't add up to one (around 0.9994). I suspect it has to do with something that happens while reading
the training file. I have yet to find the cause for this. 
