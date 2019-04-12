// Trie.hpp
// Einfache Trie-Klasse auf Vektorbasis in C++ 
// TH, 28.1.19
// Bearbeitet 26.3.19 Santichai Pornavalai

// V2: zus�tzlich: Interne Datenstrukturen, private Funktionen und Move-to-front
// V3: Key added to end state. Modified to be used as both string storing and indexing
// as well as suffix handling
// Compiler: g++

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <cstdlib> // f�r abort()
#include <vector> 


/// Maximalanzahl der Zust�nde im Trie
#define MAXSTATES   1000000

/// Trie implementiert eine Trie-Struktur auf der Basis zweier (statischer) Vektoren
class Trie
{
public:
  /// Zustandstyp des Tries
  typedef int State;
  
public:
  /// Erstellt einen leeren Trie mit Startknoten
  Trie() 
  {
    // Startzustand erzeugen (erster benutzerdef. Zustand erh�lt Index 1)
    next_free_state_index = 1;
    num_comparisions = 0;
    // Startzustand erzeugen    
    //new_state();
  }
  
  /** 
    \brief Erstellt einen Trie aus dem Inhalt der Datei 'filename' 
    \param filename Textdatei (ein Token pro Zeile, ISO-8859-kodiert)
  */
 
  Trie(std::string filename) : Trie()
  {
    // NB. Konstruktoren, die scheitern k�nnen und Fehlermeldungen ausgeben, sind generell
    // keine gute Idee. Im n�chsten Kurs verbessern wir dies.
    std::ifstream words_in(filename.c_str());
    if (words_in) {
      std::string word;
      while (std::getline(words_in,word)) {
        if (!word.empty())
          count(word);
      } // while
    }
    else std::cerr << "ERROR: Unable to open text file '" << filename << "'\n";
  }
  
  /// Destruktor
  ~Trie() 
  {
    std::cerr << "#comparisions: " << num_comparisions << "\n";
  }
  
  /** 
    \brief F�gt ein Token in den Trie ein. 
    \param token Das Token
    \return Gibt die Frequenz von 'token' zur�ck (inkl. des gerade gez�hlten)
  */
  unsigned count(std::string token)
  {
    State current_state = start_state();
    for (int i = 0; i < token.length(); ++i) {
      auto next = delta(current_state,token[i]);
      if (next == -1) {
        // Kein Folgezustand => neuen �bergang und Zielzustand erzeugen
        next = add_transition(current_state,token[i]);
        if (next == -1) {
          std::cerr << "FATAL: No room left for new transition\n";
          std::cerr << current_state;
          abort();
        }
      }
      current_state = next;
    } // for
    //current_state.
    return  increase_counter(current_state);
   // return assign_key(current_state, key);

  }


  /// adds a key to the final state
  /// used for accessing word index positions in 
  /// emmision arrays etc.
  int add_key(std::string token, int key)
  {
    State current_state = start_state();
    for (int i = 0; i < token.length(); ++i) {
      auto next = delta(current_state,token[i]);
      if (next == -1) {
        next = add_transition(current_state,token[i]);
        if (next == -1) {
          std::cerr << "FATAL: No room left for new transition\n";
          std::cerr << current_state<<std::endl;
          std::cerr <<num_types();

          abort();
        }
      }
      current_state = next;
    } // for
    //current_state.
    increase_counter(current_state);
    return assign_key(current_state, key);

  }

 
  
  


  /// same as add_key but backwards
  unsigned add_back_key(std::string token, int key){
    
    State current_state = start_state();
    for (int i =  token.length() -1; i > -1 ; --i) {
      auto next = delta(current_state,token[i]);
      if (next == -1) {
        next = add_transition(current_state,token[i]);

        if (next == -1) {
          std::cerr << "FATAL: No room left for new transition\n";
          std::cerr << current_state<<std::endl;
          std::cerr <<num_types();

          abort();
        }
      }
      current_state = next;
    } // for
    //current_state.

    increase_counter(current_state);
    return assign_key(current_state, key);

  }


  /// same as add_back_key but modified 
  unsigned add_suffix(std::string token, long long int &key){
    
    State current_state = start_state();

    for (int i =  token.length() -1; i > -1 ; --i) {
      auto next = delta(current_state,token[i]);
      if (next == -1) {
        ++key;
        next = add_transition(current_state,token[i]);
        assign_key(next, key);
        

        if (next == -1) {
          std::cerr << "FATAL: No room left for new transition\n";
          std::cerr << current_state<<std::endl;
          std::cerr <<num_types();
          abort();
        }

      }
      current_state = next;
      
    } // for
    
    assign_key(current_state, key);
    return key;
  }

/// adds suffix and all indices along the way
/// returns indices of states visited at the end
std::vector<unsigned int> add_get_suff_idx(std::string token, long long int &key){

    State current_state = start_state();
    std::vector<unsigned int> results; 

    for (int i =  token.length() -1; i > -1 ; --i) {
      auto next = delta(current_state,token[i]);
      

      if (next == -1) {
        ++key;
        next = add_transition(current_state,token[i]);
        assign_key(next, key);
        //results.push_back(states[next].key);
        

        if (next == -1) {
          std::cerr << "FATAL: No room left for new transition\n";
          std::cerr << current_state<<std::endl;
          std::cerr <<num_types();
          abort();
        }

      }
      results.push_back(states[next].key);
      current_state = next;
      
    } // for

    //increase_counter(current_state);
    
    assign_key(current_state, key);
    //results.push_back(states[current_state].key);
    return results;

}

/// gets indices of sub suffixes of a suffix
std::vector<unsigned> get_suffix_idx(std::string token){
    std::vector<unsigned> results;
    //results.reserve(5); 
    State current_state = start_state();

    for (int i =  token.length() -1; i > -1 ; --i) {
      State next = delta(current_state,token[i]);
      if (next == -1) {
        return results;
      }
      
      results.push_back(states[next].key);
      current_state = next;
    }
    return results;
}

  /// get the value from a suffix
  int get_key_suffix(std::string token ){
    State current_state = start_state();
    for (int i =  token.length() -1; i > -1 ; --i) {
      State next = delta(current_state,token[i]);
      if (next == -1) {
        return -1;
      }
      current_state = next;
    } // for
    return states[current_state].key;
  }



  int get_key(std::string token){
    State current_state = start_state();
    for (int i = 0; i < token.length(); ++i) {
      State next = delta(current_state,token[i]);
      if (next == -1) {
        return -1;
      }
      current_state = next;
    } // for
    return states[current_state].key;
  }

  /// Gibt true zur�ck, gdw. 'token' im Trie vorhanden ist
  bool in(std::string token)
  {
    return count_of(token) > 0;
  }

  int get_index(std::string token){
    return 0;
  }
  
  /// Gibt true zur�ck, gdw. 'token' im Trie vorhanden ist
  unsigned int count_of(std::string token)
  {
    State current_state = start_state();
    for (int i = 0; i < token.length(); ++i) {
      State next = delta(current_state,token[i]);
      if (next == -1) {
        return 0;
      }
      current_state = next;
    } // for
    return count_for(current_state);
  }

  /// Gibt die Anzahl der Types im Trie zur�ck
  unsigned int num_types() 
  {
    unsigned nt = 0;
    // �ber alle Zust�nde iterieren
    for (auto i = 0; i < next_free_state_index; ++i) {
      if (states[i].freq > 0)
        ++nt; 
    }
    return nt;
  }

  /// Gibt die Anzahl der Tokens im Trie zur�ck
  unsigned int num_tokens() 
  {
    unsigned nt = 0;
    // �ber alle Zust�nde iterieren
    for (auto i = 0; i < next_free_state_index; ++i) {
      nt += states[i].freq;
    }
    return nt;
  }

  /// Gibt die Anzahl der Knoten im Trie zur�ck
  unsigned int num_nodes()
  {
    return next_free_state_index;
  }

  /// Erstellt eine dot-Repr�sentation des Tries
  void as_dot(std::string dot_filename)
  {
  }
  
  /// Erstellt eine JSON-Repr�sentation des Tries
  void as_json(std::string json_filename)
  {
  }
  
  /// Erstellt eine XML-Repr�sentation des Tries
  void as_xml(std::string xml_filename)
  {
  }
  
  /// Gibt alle Elemente des Tries aus
  void print()
  {
    print_types_with_prefix("");
  }

  /// Gibt alle Types aus dem Trie aus, die mit 'pref' beginnen
  void print_types_with_prefix(std::string pref)
  {
    
  }
  
  /// Liest den Trie aus einer Bin�rdatei
  /// => C++ II
  bool read(/*...*/)
  {
    return false;
  }

  /// Schreibt den Trie im Bin�rformat in eine Datei
  /// => C++ II
  bool write(/*...*/)
  {
    return false;
  }

private: // Typen
  /// Typ eines internen Zustands
  struct InternalState {
    InternalState(unsigned int f=0, int i=-1, int k = -1) { freq = f; tr_idx = i; key = k;}
    unsigned int  freq;     ///< Assoziierte Type-Frequenz
    int           tr_idx; ///< Anfang der Adjadenzliste f�r Zustand im Vektor transitions.
    int           key; 
  };

  /// Typ eines internen �bergangs
  struct Transition {
    Transition(char s='\0', State n=-1, int nt=-1) 
    { symbol=s; next_state=n; next_tr=nt; }
    char  symbol;       ///< Symbol am �bergang
    State next_state;   ///< Folgezustand des �bergangs
    int   next_tr;      ///< N�chster Eintrag in der Adjazenzliste
  }; // Transition

private: // Funktionen
  /// Gibt den Startzustand zur�ck (momentan immer 0)
  State start_state() 
  {
    return 0;
  }
  
  /// R�ckgabe == true gdw. q ein Endzustand ist
  bool is_final(State q) 
  {
    assert(q >= 0 && q < MAXSTATES);
    return states[q].freq > 0;
  }
  
  /// Gibt den Nachfolgezustand (Determinismus!) von q mit a zur�ck, andernfalls
  /// -1 wenn dieser nicht existiert. Implementiert die Move-to-front-Strategie, d.h.
  /// gefundene �berg�nge werden an die erste Position der Adjazenzliste gesetzt.
  /// Auf diese Weise befinden sich h�ufig benutzte �berg�nge zu Beginn der Liste

  State delta(State q, char a) 
  {
    assert(q >= 0 && q < MAXSTATES);
    int last_idx = -1;
    // Lineare Suche �ber die Adjazenzliste
    for (auto idx = states[q].tr_idx; idx != -1; idx = transitions[idx].next_tr) {
      ++num_comparisions;
      // Symbol pr�fen
      if (transitions[idx].symbol == a) {
        // Ja, es ist das gesuchte => Move-front-Strategie umsetzen und Folgezustand zur�ckgeben
        if (last_idx != -1) {
          // Gefundenen �bergang "ausketten"
          transitions[last_idx].next_tr = transitions[idx].next_tr;
          // Neuer erster �bergang verweist auf alten ersten
          transitions[idx].next_tr = states[q].tr_idx;
          // Gefundenen �bergang zum ersten der Adjazenzliste machen
          states[q].tr_idx = idx;
        }
        return transitions[idx].next_state;
      }
      last_idx = idx;
    }
    // Symbol nicht gefunden
    return State(-1);
  }
  
  /// F�gt einen �bergang von Zustand q mit Symbol a zu einem NEUEN Zustand q' ein;
  /// Dieser ist auch der R�ckgabewert.
  State add_transition(State q, char a) 
  {
    if (next_free_state_index == MAXSTATES) {
      return State(-1);
    }
    assert(q >= 0 && q < MAXSTATES);
    // Neuen �bergang am Anfang der Adjazenzliste konstruieren; 
    // Liste geht am alten Beginn weiter
    transitions[next_free_state_index-1] = Transition(a,next_free_state_index,states[q].tr_idx);
    // Im Zustandsvektor an der Stelle q Beginn der neuen Liste festlegen
    states[q].tr_idx = next_free_state_index-1;
    // Zustandsindex zur�ckgeben und "anschlie�end" hochz�hlen
    return next_free_state_index++;
  }
  
  /// Erh�ht den Token-Z�hler f�r Zustand q. Gibt den erh�hten Wert zur�ck
  unsigned int increase_counter(State q) 
  {
    assert(q >= 0 && q < MAXSTATES);
    return ++states[q].freq;
  }

    unsigned int assign_key(State q, int key) 
  {
    assert(q >= 0 && q < MAXSTATES);
    return states[q].key = key;
  }
  /// Gibt den Token-Z�hler f�r Zustand q zur�ck
  unsigned int count_for(State q) 
  {
    assert(q >= 0 && q < MAXSTATES);
    return states[q].freq;
  }

private: // Instanzvariablen
  InternalState states[MAXSTATES];
  Transition    transitions[MAXSTATES-1];
  unsigned int  next_free_state_index;
  unsigned int  num_comparisions;
}; // Trie
