// driver file for HMM POS tagger
// Santi(chai) Pornavalai
// 1.4.19
// tested on g++ -O2 
#include <iostream> 
#include "hmm.hpp"
#include "Trie.hpp"
#include <chrono> 

int main(int argc, char** argv){
  HMM* hmm = new HMM(); 

  if (argc < 3 or argc >4){ 
    std::cerr << "ERROR Invalid Number of Arguments"<<std::endl;
    std::cerr << argv[0]<<" <TRAIN FILE> <TEST FILE>"<<std::endl;
    exit(2);}
  
  auto start = std::chrono::high_resolution_clock::now();
  std::cerr <<"training"<<std::endl;
  hmm->read_text(argv[1]);
  hmm->normalize();
  std::cerr<<"tagging"<<std::endl;
  hmm->tag_all(argv[2]);
  auto finish = std::chrono::high_resolution_clock::now();

  
  std::chrono::duration<double> elapsed = finish - start;
  std::cerr << "Elapsed time: " << elapsed.count() << " s \n";

  return 0;
}