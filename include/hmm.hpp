/*
HMM class for POS tagging
Santi(chai)	Pornavalai
30.3.19
*/

#include<iostream> 
#include<string>
#include<vector> 
#include<cassert> 
#include<fstream>  

#include <functional> // for accumulate. makes code more readable
#include <numeric>

#include "Trie.hpp"
#include "word_map.hpp"
#include <cstdlib> // for abort
#include <cmath> //floating point comparison
#include <map> // for suffix debugging... promise

#pragma once
//#define NDEBUG


///HMM Class for training, tagging and suffix handling
class HMM{


public:

	/// default constructor
	/// builds empty matrices and initializes index Tries
	HMM(){


		word2id = new Trie(); 
		suffix2id = new Trie();
		tag2id = new Trie();
		id2word = new WordMap();
		id2tag = new WordMap();	
		initialize_weights();

		};


	/// TODO:
	 HMM(const HMM &hmm) 
   { 
      std::cout<<"Copy constructor not yet implemented"<<std::endl; 
   };

   HMM& operator = (const HMM &t) 
   { 
      std::cout<<"Assignment constructor not yet implemented"<<std::endl; 
   } 

	~HMM(){
		delete word2id; 
		delete suffix2id;
		delete tag2id;
		delete id2word;
		delete id2tag;
		};

	

	public: 

	/// reads file
	/// \param fname  name of file to read and tag
	void tag_all(std::string fname){
			
			std::vector<std::vector<std::string>> in_vec =  read_sents(fname);
			for(unsigned int i = 0; i < in_vec.size(); ++i){
				decode(in_vec[i]);
				std::cout<<std::endl;
			}
	}



	/// takes vector of strings and tags it using viterbi algorithm
	/// unknown words are handled using a separate probability for suffixes
	/// \param in_string  input vector of strings
	std::vector<std::string> decode(std::vector<std::string> in_string){

		unsigned int in_size = in_string.size();
		double word_prob; 

		DoubleMatrix forward_table(in_size, std::vector<double>(next_free_tag,0.0));
		IntMatrix bp_table(in_size, std::vector<int>(next_free_tag,0));

		std::string first_word = in_string[0];
		int first_id = word2id->get_key(first_word);
		
		// initialize viterbi trellis
		for(unsigned i = 0; i < next_free_tag; i++){

			 forward_table[0][i] = initial_prob[i]*emmision_prob_weighted(first_word,i);
		}

		//viterbi algorithm

    // for word
		for(unsigned k = 1; k < in_size; k++){           
        // for tag
				for(unsigned i = 0; i < next_free_tag; i++){
					word_prob = emmision_prob_weighted(in_string[k], i);
					// for prev tag
					for(unsigned j= 0; j < next_free_tag; j++){
						double trans = tag_prob_id(i,j) ;
						//find arg max
						if (trans * word_prob* forward_table[k-1][j] > forward_table[k][i] ){                   
						forward_table[k][i] = trans * word_prob * forward_table[k-1][j];
						bp_table[k][i] = j;
                        
								}
							}
						}
				}

			// backtrack 
			std::vector<std::string> results;
			results = trace_back(bp_table, forward_table);
			for(int i = 0; i  <in_size; i++){
			std::cout<< in_string[i]<<" " <<results[i] <<std::endl;
			}

		return results;
	}



	
/// reads sentences from file name and stores them in a matrix of strings
/// \param fname name of file
/// \ param delim delimiter used. default is space
/// \returns vector of string vectors
std::vector<std::vector<std::string>> read_sents(std::string fname, std::string delim = " "){

	std::cerr<<"opened file for reading sents" <<std::endl;
  std::ifstream f(fname);

	if(!f.good()){
		std::cerr << "File does not exis/or is corrupted."<<std::endl;
		exit(2);
	}

	std::string word="";
	std::string tag="";
	std::string tag_1="";
	std::string tag_2=""; 
	std::string line; 

	std::vector <std::vector <std::string>>  results; 
	std::vector<std::string> curr_sent;

	while(getline(f,line)){
		//std::cout<<line; 

		if(line.length()<1){
			results.push_back(curr_sent);
			curr_sent = {};
		}


		else{
			unsigned split = line.find(delim);
			word = line.substr(0,split);
			tag = line.substr(split+1,line.length());
			curr_sent.push_back(word);
			}
		}


		results.push_back(curr_sent);
		return results;
	}



	/// initialize matrices used
	/// size is determined by private variables
	/// it is ok to have smaller matrices than needed
	/// they grow when about to run out of space
	/// also attempts to catch bad alloc
	void initialize_weights(){
		// try/catch clause for bad_alloc when given size is too large

		try{

		bigram_prob = std::vector<std::vector<double>>
		(tag_size, std::vector<double>(tag_size,0.0));

		unigram_prob = std::vector<double>(tag_size,0.0);
		initial_prob = std::vector<double>(tag_size,0.0);
		emmision_prob = std::vector<std::vector<double>>
		 (tag_size, std::vector<double>(word_size,0.0));

		suffix_prob = std::vector<std::vector<double>>
		(tag_size, std::vector<double>(word_size,0.0));

		}

		catch(std::bad_alloc & ba){
			std::cerr << "bad alloc caught: " <<ba.what()<<"\n";
			std::cerr << "aborting HMM"<<"\n";
			exit(3);
		}

		
	}


	/// read text and fills up emmision, transition and suffix probability vectors
	/// \param fname file name in <word>   <tag> form. separated by exactly 3 spaces... sorry
	/// \param delim  delimiter . default to 3 spaces
	bool read_text(std::string fname, std::string delim = "   "){

	std::cerr<<"opened file" <<std::endl;
  std::ifstream f(fname);

	if(!f.good()){
		std::cerr << "File does not exis/or is corrupted."<<std::endl;
		exit(2);
	}
	// default constructor of std::string is empty string
	std::string word, tag,  line; 
	int idx_prev_tag = 0;
	bool is_start = false;


	while(getline(f,line)){

		if(line.length()<1){
			is_start = true;
		}

		else{

			unsigned split = line.find(delim);
			
			word = line.substr(0,split);
			tag = line.substr(split+1,line.length());
			//std::cout<< word <<" "<<tag<<std::endl;

			int idx_word = word_in_idx(word);
			int idx_tag = tag_in_idx(tag);
			add_tag_seq(idx_tag,idx_prev_tag);
			add_word_tag(idx_word,idx_tag);

			std::vector<unsigned int> idxes_suffix = add_suffix(word);
			increment_suffix(idxes_suffix,idx_tag);
			//std::cout<< idx_tag<<"   " << idx_prev_tag << std::endl;

			if(is_start){
				initial_prob[idx_tag] += 1;
				is_start = false;
			}
			
			idx_prev_tag = idx_tag;

			}
		}
		f.close();
		return true;
	}


/// used to debug hmm weights
/// epsilon determines the precision of the comparison
void verify(double epsilon = 0.000001){
	double sum = 0;
	for(unsigned int i = 0; i < next_free_tag; ++i){
		sum = 0;
		double sum_uni = 0;
		double sum_emi = 0;
		double sum_suffix = 0;
		

		for(unsigned int j = 0; j < next_free_tag; ++j){

				sum += bigram_prob[i][j];
				sum_uni += unigram_prob[j];
			}

			if(not essentiallyEqual(sum,1, epsilon)){
				std::cerr<< id2tag ->get_word(i)<<std::endl;
				std::cerr<< "in bigram ";
				std::cerr<<sum<<std::endl;
			}

			if(not essentiallyEqual(sum_uni,1, epsilon)){
				std::cerr<< id2tag ->get_word(i)<<std::endl;
				std::cerr<<sum_uni<<std::endl;
			}

			for(unsigned int j = 0; j < next_free_word; ++j){
				sum_emi += emmision_prob[i][j];
			}

			if(not essentiallyEqual(sum_emi,1, epsilon)){
				std::cerr<< id2tag->get_word(i)<<std::endl;
				std::cerr<<sum_emi<<std::endl;
			}

			for(unsigned int j = 0; j < next_free_suffix; ++j){
				sum_suffix += suffix_prob[i][j];
			}

			if(not (essentiallyEqual(sum_suffix,1, epsilon) or essentiallyEqual(sum_suffix,0, epsilon) ) ){
				std::cerr<< id2tag->get_word(i)<<std::endl;
				std::cerr<< "in suffix ";

				std::cerr<<sum_suffix<<std::endl;
			}


		}
	}



/// normalize all the probability matrices
void normalize(){

	double sum;
	double sum_emi;
	double sum_suffix;

	double sum_unigrams = std::accumulate(unigram_prob.begin(), unigram_prob.end(), 0.0);
	double sum_initial = std::accumulate(initial_prob.begin(), initial_prob.end(), 0.0);

	

	for(unsigned int i = 0; i < next_free_tag; ++i){
		sum = unigram_prob[i];
		sum_emi = std::accumulate(emmision_prob[i].begin(), emmision_prob[i].end(), 0.0);
		sum_suffix = std::accumulate(suffix_prob[i].begin(), suffix_prob[i].end(), 0.0);

		//tag & initial
		unigram_prob[i] /= sum_unigrams;
		initial_prob[i] /= sum_initial;

		// tag-tag
		for(unsigned int j = 0; j < next_free_tag; ++j){
			bigram_prob[i][j] /= sum;
			}
		// tag - word
		for(unsigned int j = 0; j < next_free_word; ++j){
			emmision_prob[i][j] = (emmision_prob[i][j] + smoothing_k) / (sum_emi + smoothing_k);
			}
		// tag - suffix
		if (sum_suffix>0.0){
		for(unsigned int j = 0; j < next_free_suffix; ++j){
			suffix_prob[i][j] /=  sum_suffix;
			}
		
		}
		}
	}



	/// prints out emmisions, transition, suffix, unigram etc. to standard output
void print_weights(){
	std::cout<< "# TRANSITION WEIGHTS";
	for(int i = 0; i < next_free_tag; ++i){
		for(int j = 0; j < next_free_tag; ++j){
			double weight = bigram_prob[i][j];

			if ((weight > 0)) {
				std::cout<< id2tag->get_word(i)<< "\t"<< id2tag->get_word(j)<<"\t"
				<<weight<<std::endl;
			}
						
		}
	}

	std::cout<< "# UNIGRAM WEIGHTS";
	for(int i = 0; i < next_free_tag; ++i){
		
			double weight = unigram_prob[i];

			if ((weight > 0)) {
				std::cout<< id2tag->get_word(i) << "\t"
				<<weight<<std::endl;
			}
							
	}

		std::cout<< "# INITIAL WEIGHTS";
	for(int i = 0; i < next_free_tag; ++i){
		
			double weight = initial_prob[i];

			if ((weight > 0)) {
				std::cout<< id2tag->get_word(i) << "\t"
				<<weight<<std::endl;
			}
							
	}


	std::cout<< "# EMMISION WEIGHTS";
	for(int i = 0; i < next_free_tag; ++i){
		for(int j = 0; j <= next_free_word; ++j){
			double weight = emmision_prob[i][j];

			 if(weight > smoothing_k){
				std::cout<< id2word ->get_word(j) << "\t"<< id2tag->get_word(i)<<"\t"
				<<weight<<std::endl;
			 }
						
		}
	}

	std::cout<< "# Suffix WEIGHTS";
	for(int i = 0; i < next_free_tag; ++i){
		for(int j = 0; j <= next_free_suffix; ++j){
			double weight = suffix_prob[i][j];

			 if(weight > 0.0){
				std::cout<< id2suffix[j] << "\t"<<id2tag->get_word(i)<<"\t"
				<<weight<<std::endl;
			 }
						
		}
	}

}

private:

// typedef doesn't work here. sorry.
/// extract path from viterbi trellis and backpointer table
inline	std::vector<std::string> trace_back(std::vector<std::vector<int>> &bp_table,std::vector<std::vector<double>> &forward_table){

	double max = 0; 
	int best_pos = 0; 
	std::vector<std::string> results;
	unsigned in_size = bp_table.size();
	results.resize(in_size);


	//get max
	for(unsigned i = 0; i < next_free_tag; i++){
			if(forward_table[in_size-1][i] > max) {
				max =forward_table[in_size-1][i];
				best_pos = i;  
			}
		}	

	// trace_back
		for(int i = in_size-1; i  != -1; i--){
			 int bp = bp_table[i][best_pos];
				results[i]= id2tag->get_word(best_pos);

			best_pos = bp;	
		}

	return results;

	}


	
/// wrapper to get index of a word
/// if nout found adds word to trie
inline int word_in_idx(std::string word){
	int idx = word2id -> get_key(word);
	if(idx < 0){
			//add word to id maps
			word2id -> add_key(word, next_free_word);
			id2word -> add_word(word);
			++next_free_word;
			return next_free_word ;
			}
	return idx;

}
/// wrapper for tag getter
/// if tag not found adds tag
/// if found returns index
inline int tag_in_idx(std::string tag){
		if(tag2id -> get_key(tag) < 0){
		tag2id -> add_key(tag, next_free_tag);
	//	id2tag.insert(std::pair< int, std::string>( next_free_tag,tag));

		id2tag ->add_word(tag);
		++next_free_tag;
		// work around. 
		return next_free_tag - 1;
		}
	return tag2id ->get_key(tag);
}



/// gets emmision weight if found
/// if not switches to suffix handling
inline double emmision_prob_weighted ( std::string word, unsigned tag_id){
	// temp fix
	double word_prob;
	int word_idx = word2id -> get_key(word);

	
	if(word_idx < 0){
		word_prob = get_suffix_prob(get_suffix(word),tag_id);

	}
	else{
		word_prob = emmision_prob[tag_id][word_idx];

	}
	

	return word_prob;
}

/// gets the interpolated probability of tag uni/bigram
inline double tag_prob(std::string tag, std::string prev_tag,double lambda_bi = 0.95){


	int idx_tag = tag2id ->get_key(tag);
	int idx_prev = tag2id ->get_key(prev_tag);
	double lambda_uni = 1 - lambda_bi;
	return lambda_bi* bigram_prob[idx_tag][idx_prev] + 
			lambda_uni*unigram_prob[idx_prev];
}


/// same as tag prob but just uses id instead of string
inline double tag_prob_id(unsigned idx_tag, unsigned idx_prev, double lambda_bi = 1){

	double lambda_uni = 1 - lambda_bi;
	return lambda_bi* bigram_prob[idx_tag][idx_prev] + 
			lambda_uni*unigram_prob[idx_prev];
}

/// adds tag sequence to transition probabilities
/// resizes if limit is reached
inline bool add_tag_seq(int tag, int prev_tag){
	// bigram[crr][prev]

	//assert(prev_tag > tag_size -1);
	if (tag > tag_size-1 ){
		std::cerr<<"overflowing";

		++tag_size;
		std::cerr<<tag_size;
		std::vector<double> tmp_tag_vec = std::vector<double>(tag_size,0.0);
		
		tmp_tag_vec[prev_tag] += 1; 
		bigram_prob.push_back(tmp_tag_vec);

		std::vector<double> tmp_word_vec = std::vector<double>(word_size,0.0);
		emmision_prob.push_back(tmp_word_vec);

		std::vector<double> tmp_suffix_vec = std::vector<double>(word_size,0.0);
		suffix_prob.push_back(tmp_suffix_vec);

		for(int i = 0; i < bigram_prob.size(); ++i){
			bigram_prob[i].push_back(0.0);
		}
		
		unigram_prob.push_back(0.0);
		initial_prob.push_back(0.0);

		return false;
	}

	bigram_prob[tag][prev_tag] += 1.0;
	unigram_prob[prev_tag]+= 1.0;
	return true;

	}


/// adds suffix to trie
/// returns a vector of indices indicating where to find each of the sub suffixes
/// in the suffixes in the trie...
inline std::vector<unsigned int> add_suffix(std::string word ){

	std::vector<unsigned int> results; 

	if (word.length() > 8){
				word = word.substr(word.length()-suff_cutoff,word.length());
			
	results =  suffix2id -> add_get_suff_idx(word,next_free_suffix);
	for(int i = 0; i < results.size(); ++i){

		if(id2suffix.find(results[i]) == id2suffix.end()){
		id2suffix.insert(std::pair< int, std::string>( results[i],word.substr(word.length()-i-1,word.length())));

				}

		}
	}
	
	return results;
}


/// gets suffix but does so without adding
inline std::vector<unsigned int> get_suffix(std::string word ){
	std::vector<unsigned int> results; 

	if (word.length() > suff_cutoff){
				word = word.substr(word.length()-suff_cutoff,word.length());
			}
	
	results =  suffix2id -> get_suffix_idx(word);
	return results;
}

/// increment the suffix in the suffix probability matrix
inline void increment_suffix(std::vector<unsigned int> indices , unsigned tag_id){
	for(int i =0; i<indices.size(); i++){
		if(indices[i] > word_size){
				for(int j = 0; j < suffix_prob.size(); ++j){
				suffix_prob[i].push_back(0.0);
			}

		}
		suffix_prob[tag_id][indices[i]] += suff_weights[i];

	}

}


/// gets the sum off all the the sub-suffixes probabilities given a tag
inline double get_suffix_prob(std::vector<unsigned int> indices , unsigned tag_id){
	double accu = 0;
	for(int i =0; i<indices.size(); i++){
		accu += suffix_prob[tag_id][indices[i]];
	}
	return accu; 
}

inline bool add_word_tag(int word_pos, int tag_pos){
	if (word_pos > word_size -1){
		++word_size;
				std::abort();

		for(int i = 0; i < emmision_prob.size(); ++i){
			emmision_prob[i].push_back(0.0);
		}

	}
	emmision_prob[tag_pos][word_pos] += 1;
	return true;
}





/// floating point comparison coutersy of D.Knuth.
bool essentiallyEqual(double a, double b, double epsilon)
{
    return std::fabs(a - b) <= ( (std::fabs(a) > std::fabs(b) ? std::fabs(b) : std::fabs(a)) * epsilon);
}

private:

/// predicted number of words
/// a good estimation makes training a bit quicker
int word_size =100000;

/// predicted number of tags
int tag_size = 60;

/// unigram probability
std::vector<double> unigram_prob;

/// transition probabilities
std::vector<std::vector<double>> bigram_prob;

///start probs
std::vector<double> initial_prob;

/// emmision probabilities
std::vector<std::vector<double>> emmision_prob;
/// suffix probabilities
std::vector<std::vector<double>> suffix_prob;

// index getters
Trie* word2id  ;
Trie* suffix2id ;
Trie* tag2id; 

// inverse of the above
WordMap* id2word = new WordMap();
WordMap* id2tag = new WordMap();
std::map<int, std::string> id2suffix;


// max indices for keeping track
unsigned next_free_word  = 0;
unsigned next_free_tag = 0;
long long int next_free_suffix = 0; 

// obsolete after suffix handling
double smoothing_k = 0;
// other parameters for suffix handling
int suff_cutoff = 5; 
std::vector<double> suff_weights = {0.5,0.5,1,1,1,1};
// some type defs
typedef std::vector<std::vector<double>>  DoubleMatrix;
typedef std::vector<std::vector<int>>  IntMatrix;
typedef std::vector<std::vector<std::string>>  StringMatrix;

};