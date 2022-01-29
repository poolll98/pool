#include "split_and_merge_algorithm.h"
#include <iostream>
#include <algorithm>
#include <random>

void SplitAndMergeAlgorithm::print_startup_message() const {
  std::string msg = "Running Split and Merge algorithm with " +
                    bayesmix::HierarchyId_Name(unique_values[0]->get_id()) +
                    " hierarchies, " +
                    bayesmix::MixingId_Name(mixing->get_id()) + " mixing...";
  std::cout << msg << std::endl;
}



void SplitAndMergeAlgorithm::compute_S(const unsigned int i, const unsigned int j) {
    unsigned int lengthS;
    for(int k; k<allocations.size(); k++){
        if(allocations[k]==allocations[i]||allocations[k]==allocations[j])
            lengthS++;
    }
    lengthS=lengthS-2;
    S={}; //Necessary, since it needs to be zero at each iter, not only at the beginning
    S=S.resize(lengthS,0);
    unsigned int index=0;
    for (int k = 0; k < allocations.size(); ++k) {
        if ((allocations[k]==allocations[i]||allocations[k]==allocations[j])&& (k!=i) && (k!=j)) {
            if (index>=lengthS)
                std::cerr<< "Index out of bounds. index="<<index
                    <<",lengthS="<<lengthS<< " allocations="<<allocations<<", i="<<i<<", j="<<j<<std::endl;
            S[index]=k;
            index++;
        }
    }
}



std::vector<unsigned int> SplitAndMergeAlgorithm::compute_C_launch(const unsigned int i, const unsigned int j){
    if (allocations[i]==allocations[j]) {
        LabI = *max_element(allocations.begin(), allocations.end()) + 1;
    }else{
        LabI=allocations[i];
    }
    std::vector<unsigned int>cl(S.size(),LabI);
    std::default_random_engine generator;
    std::bernoulli_distribution bdistr(0.5);
    for (int k = 0; k < S.size(); ++k) {
        if (bdistr(generator))
            cl[k]=allocations[j];
    }
    return cl;
}



void SplitAndMergeAlgorithm::split_or_merge(std::vector<unsigned int>& cl, const unsigned int i, const unsigned int j){
    if(allocations[i]==allocations[j]) { 
      LabI=*(std::max_element(allocations.begin(), allocations.end()));
      std::vector<unsigned int> clSplit (allocations.size()); #we could initialize the vector to LAbI
      clSplit[i]=LabI;
      clSplit[j]=allocations[j];
      unsigned int CountLabI=0;
      unsigned int CountLabJ=0;
      const double q=restricted_GS(cl,i,j,1);
      unsigned int z=0;
      for(unsigned int i=0; i < clSplit.size(); i++){
          if((z<S.size())and(i==S[z])){
            if(cl[z]==LabI) CountLabI++;
            else CountLabJ++;
            
            clSplit[i]=cl[z];
            z++;
                                      }
          else{
            clSPlit[i]=allocations[i];
              }
                                            
                                            }
      const double p1=1/q;
      const double p2=factorial(CountLabI-1)*factorial(CountLabJ-1)/(S.size()+2-1)*hierarchy.alpha;
      
      
                                                                                              }
bool accepted_proposal(const double acRa) const{
    std::default_random_engine generator;
    std::uniform_real_distribution UnifDis(0.0, 1.0);
    return (UnifDis(generator)<=acRa);
                                                }
  # standard Gibbs Sampling
void restricted_GS(std::vector<unsigned int>& cl, const unsigned int i, 
                   const unsigned int j, double &res_prod){ #è stata messa _const non so perchè 
  for(unsigned int i=0; i<S.size(); i++){
    LabI=*(std::max_element(allocations.begin(), allocations.end())); #bisogna mettere LabI come _private
    p_i = ComputeRestrGSProbabilities(cl, i, j, z, 'i');
    p_j = ComputeRestrGSProbabilities(cl, i, j, z, 'j');
    p   = p_i/(p_i+p_j);  
    cl[i]= (accepted_proposal(p)) ? LabI : false;
                                         }                                                        
                                                           }
      
      
      
      
      
      
      
      
      
      
      
