// integer to word map
// using char vectors
// Santi(chai) Pornavalai
// 20.3.19

#include <vector> 
#include <string>
#include <iostream> 
#include <stdexcept>

///  Class for mapping integers to strings in a space
/// efficient way
class WordMap{

public: 

WordMap(unsigned size_word=250000, unsigned size_char=80000){
// avg. word length for german * avg vocab size in corpus

character_array.reserve(250000);
word_ids.reserve(80000);

}

~WordMap(){

}

/// free space
void shrink(){
    character_array.shrink_to_fit();
    word_ids.shrink_to_fit();

}

/// adds word to index
void add_word(std::string word){
    word_ids.push_back(curr_char_pos);
    for(unsigned i=0; i < word.size(); i++){
        character_array.push_back(word[i]);
        ++ curr_char_pos;
    }
    character_array.push_back(0);
    ++curr_char_pos;
}

/// gets word to index
std::string get_word(unsigned index){
    std::string res; 
    if(index > word_ids.size() - 1) throw std::out_of_range("Indexerror");
    unsigned pos = word_ids[index]; 

    for( int i = pos; int(character_array[i]) != 0; ++i){
        res += character_array[i];
    }
    
    return res; 
}

/// checks if word is in index given index
bool word_in(unsigned indx){
    if (indx > curr_num){
        return false;
    }
    return true;
}
/// print out words for debugging
void print(){

    for(unsigned i=0; i < curr_char_pos; ++i){
        if (int(character_array[i]) == 0) std::cout<< "NULL" <<std::endl;
        else  std::cout<< character_array[i] << std::endl;
    }
}

private: 

std::vector<int> word_ids; 
std::vector<char> character_array; 
unsigned curr_num = 0;
unsigned curr_char_pos = 0; 

};