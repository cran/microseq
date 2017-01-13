#include <fstream>
#include <sstream>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List read_fastq(std::string path) { // for standard fastq, not Sanger
  std::vector<std::string> header, sequence, quality;
  std::ifstream t(path.c_str());
  std::stringstream ss;
  ss << t.rdbuf();
  std::string to;

  while(std::getline(ss,to,'\n')){
    header.push_back(to);
    std::getline(ss,to,'\n');
    sequence.push_back(to);
    std::getline(ss,to,'\n'); // +
    std::getline(ss,to,'\n'); // Quality
    quality.push_back(to);
  }
  // return Rcpp::DataFrame::create( Rcpp::Named("Header")= header, Rcpp::Named("Sequence") = sequence);
  return Rcpp::List::create(Rcpp::Named("Header") = header, Rcpp::Named("Sequence") = sequence, Rcpp::Named("Quality") = quality);
}

// [[Rcpp::export]]
Rcpp::List read_fastq_Sanger(std::string path) { // for fastq with Sanger compatibility
  std::vector<std::string> header, sequence, quality;
  std::ifstream t(path.c_str());
  std::stringstream ss;
  ss << t.rdbuf();
  std::string to;
  unsigned int n_lin = 1, n_seq = 0;
  
  // First round
  std::getline(ss,to,'\n');
  header.push_back(to);
  std::getline(ss,to,'\n');
  sequence.push_back(to);
  while(std::getline(ss,to,'\n')){
    if(to == "+" || (to.substr(0,1) == "+" && to.substr(1, std::string::npos) == header[n_seq].substr(1,std::string::npos))){
      break;
    } else { // Count number of lines used for current sequence
      sequence[n_seq].append(to);
      n_lin++;
    }
  }

  while(std::getline(ss,to,'\n')){
    // Read number of lines used for current sequence (quality)
    quality.push_back(to);
    for(unsigned int i = 1; i<n_lin; ++i){
      std::getline(ss,to,'\n'); 
      quality[n_seq].append(to);
    }
    n_lin = 1;

    // Store header and sequence
    std::getline(ss,to,'\n');
    if(to == ""){
      break;
    } else {
      n_seq++;
    }
    header.push_back(to);
    std::getline(ss,to,'\n');
    sequence.push_back(to);
    while(std::getline(ss,to,'\n')){
      if(to == "+" || (to.substr(0,1) == "+" && to.substr(1, std::string::npos) == header[n_seq].substr(1,std::string::npos))){
        break;
      } else { // Count number of lines used for current sequence
        sequence[n_seq].append(to);
        n_lin++;
      }
    }
  }
  // return Rcpp::DataFrame::create( Rcpp::Named("Header")= header, Rcpp::Named("Sequence") = sequence);
  return Rcpp::List::create(Rcpp::Named("Header") = header, Rcpp::Named("Sequence") = sequence, Rcpp::Named("Quality") = quality);
}

// [[Rcpp::export]]
bool write_fastq(std::vector<std::string> header, std::vector<std::string> sequence, std::vector<std::string> quality, std::string path) { // for standard FASTA format
  std::ofstream t(path.c_str());
  for(unsigned int j=0; j<header.size();++j){
    t << header[j] << std::endl;
    t << sequence[j] << std::endl;
    t << '+' << std::endl;
    t << quality[j] << std::endl;
  }
  t.close();
  return(true);
}


// [[Rcpp::export]]
Rcpp::List read_fasta(std::string path) { // for standard FASTA format
  std::vector<std::string> header, sequence;
  std::ifstream t(path.c_str());
  std::stringstream ss;
  ss << t.rdbuf();
  std::string to;
  unsigned int n_seq = 0;
  
  std::getline(ss,to,'\n');
  header.push_back(to.substr(1));
  std::getline(ss,to,'\n');
  sequence.push_back(to);
  while(std::getline(ss,to,'\n')){
    if(to.substr(0,1) == ">"){
      ++n_seq;
      header.push_back(to.substr(1));
      std::getline(ss,to,'\n');
      sequence.push_back(to);
    } else {
      sequence[n_seq].append(to);
    }
  }
  // return Rcpp::DataFrame::create( Rcpp::Named("Header")= header, Rcpp::Named("Sequence") = sequence);
  return Rcpp::List::create(Rcpp::Named("Header") = header, Rcpp::Named("Sequence") = sequence);
}

// [[Rcpp::export]]
bool write_fasta(std::vector<std::string> header, std::vector<std::string> sequence, std::string path, int width) { // for standard FASTA format
  std::ofstream t(path.c_str());
  for(unsigned int j=0; j<header.size();++j){
    t << ">" << header[j] << std::endl;
    for (unsigned i = 0; i < sequence[j].length(); i += width) {
      t << sequence[j].substr(i, width) << std::endl;
    }
  }
  t.close();
  return(true);
}
