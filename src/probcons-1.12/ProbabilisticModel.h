/////////////////////////////////////////////////////////////////
#include <cstring>
// ProbabilisticModel.h
#include <cstring>
//
#include <cstring>
// Routines for (1) posterior probability computations
#include <cstring>
//              (2) chained anchoring
#include <cstring>
//              (3) maximum weight trace alignment
#include <cstring>
/////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
#ifndef PROBABILISTICMODEL_H
#include <cstring>
#define PROBABILISTICMODEL_H
#include <cstring>

#include <cstring>
#include <list>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstring>
#include "SafeVector.h"
#include <cstring>
#include "ScoreType.h"
#include <cstring>
#include "SparseMatrix.h"
#include <cstring>
#include "MultiSequence.h"
#include <cstring>

#include <cstring>
using namespace std;
#include <cstring>

#include <cstring>
const int NumMatchStates = 1;                                    // note that in this version the number
#include <cstring>
                                                                 // of match states is fixed at 1...will
#include <cstring>
                                                                 // change in future versions
#include <cstring>
const int NumMatrixTypes = NumMatchStates + NumInsertStates * 2;
#include <cstring>

#include <cstring>
/////////////////////////////////////////////////////////////////
#include <cstring>
// ProbabilisticModel
#include <cstring>
//
#include <cstring>
// Class for storing the parameters of a probabilistic model and
#include <cstring>
// performing different computations based on those parameters.
#include <cstring>
// In particular, this class handles the computation of
#include <cstring>
// posterior probabilities that may be used in alignment.
#include <cstring>
/////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
class ProbabilisticModel {
#include <cstring>

#include <cstring>
  float initialDistribution[NumMatrixTypes];               // holds the initial probabilities for each state
#include <cstring>
  float transProb[NumMatrixTypes][NumMatrixTypes];         // holds all state-to-state transition probabilities
#include <cstring>
  float matchProb[256][256];                               // emission probabilities for match states
#include <cstring>
  float insProb[256][NumMatrixTypes];                      // emission probabilities for insert states
#include <cstring>

#include <cstring>
 public:
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ProbabilisticModel()
#include <cstring>
  //
#include <cstring>
  // Constructor.  Builds a new probabilistic model using the
#include <cstring>
  // given parameters.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  ProbabilisticModel (const VF &initDistribMat, const VF &gapOpen, const VF &gapExtend,
#include <cstring>
                      const VVF &emitPairs, const VF &emitSingle){
#include <cstring>

#include <cstring>
    // build transition matrix
#include <cstring>
    VVF transMat (NumMatrixTypes, VF (NumMatrixTypes, 0.0f));
#include <cstring>
    transMat[0][0] = 1;
#include <cstring>
    for (int i = 0; i < NumInsertStates; i++){
#include <cstring>
      transMat[0][2*i+1] = gapOpen[2*i];
#include <cstring>
      transMat[0][2*i+2] = gapOpen[2*i+1];
#include <cstring>
      transMat[0][0] -= (gapOpen[2*i] + gapOpen[2*i+1]);
#include <cstring>
      assert (transMat[0][0] > 0);
#include <cstring>
      transMat[2*i+1][2*i+1] = gapExtend[2*i];
#include <cstring>
      transMat[2*i+2][2*i+2] = gapExtend[2*i+1];
#include <cstring>
      transMat[2*i+1][2*i+2] = 0;
#include <cstring>
      transMat[2*i+2][2*i+1] = 0;
#include <cstring>
      transMat[2*i+1][0] = 1 - gapExtend[2*i];
#include <cstring>
      transMat[2*i+2][0] = 1 - gapExtend[2*i+1];
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // create initial and transition probability matrices
#include <cstring>
    for (int i = 0; i < NumMatrixTypes; i++){
#include <cstring>
      initialDistribution[i] = LOG (initDistribMat[i]);
#include <cstring>
      for (int j = 0; j < NumMatrixTypes; j++)
#include <cstring>
        transProb[i][j] = LOG (transMat[i][j]);
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // create insertion and match probability matrices
#include <cstring>
    for (int i = 0; i < 256; i++){
#include <cstring>
      for (int j = 0; j < NumMatrixTypes; j++)
#include <cstring>
        insProb[i][j] = LOG (emitSingle[i]);
#include <cstring>
      for (int j = 0; j < 256; j++)
#include <cstring>
        matchProb[i][j] = LOG (emitPairs[i][j]);
#include <cstring>
    }
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeForwardMatrix()
#include <cstring>
  //
#include <cstring>
  // Computes a set of forward probability matrices for aligning
#include <cstring>
  // seq1 and seq2.
#include <cstring>
  //
#include <cstring>
  // For efficiency reasons, a single-dimensional floating-point
#include <cstring>
  // array is used here, with the following indexing scheme:
#include <cstring>
  //
#include <cstring>
  //    forward[i + NumMatrixTypes * (j * (seq2Length+1) + k)]
#include <cstring>
  //    refers to the probability of aligning through j characters
#include <cstring>
  //    of the first sequence, k characters of the second sequence,
#include <cstring>
  //    and ending in state i.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  VF *ComputeForwardMatrix (Sequence *seq1, Sequence *seq2) const {
#include <cstring>

#include <cstring>
    assert (seq1);
#include <cstring>
    assert (seq2);
#include <cstring>

#include <cstring>
    const int seq1Length = seq1->GetLength();
#include <cstring>
    const int seq2Length = seq2->GetLength();
#include <cstring>

#include <cstring>
    // retrieve the points to the beginning of each sequence
#include <cstring>
    SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
#include <cstring>
    SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
#include <cstring>

#include <cstring>
    // create matrix
#include <cstring>
    VF *forwardPtr = new VF (NumMatrixTypes * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
#include <cstring>
    assert (forwardPtr);
#include <cstring>
    VF &forward = *forwardPtr;
#include <cstring>

#include <cstring>
    // initialization condition
#include <cstring>
    forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] = 
#include <cstring>
      initialDistribution[0] + matchProb[(unsigned char) iter1[1]][(unsigned char) iter2[1]];
#include <cstring>
   
#include <cstring>
    for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
      forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] = 
#include <cstring>
	initialDistribution[2*k+1] + insProb[(unsigned char) iter1[1]][k];
#include <cstring>
      forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] = 
#include <cstring>
	initialDistribution[2*k+2] + insProb[(unsigned char) iter2[1]][k]; 
#include <cstring>
    }
#include <cstring>
    
#include <cstring>
    // remember offset for each index combination
#include <cstring>
    int ij = 0;
#include <cstring>
    int i1j = -seq2Length - 1;
#include <cstring>
    int ij1 = -1;
#include <cstring>
    int i1j1 = -seq2Length - 2;
#include <cstring>

#include <cstring>
    ij *= NumMatrixTypes;
#include <cstring>
    i1j *= NumMatrixTypes;
#include <cstring>
    ij1 *= NumMatrixTypes;
#include <cstring>
    i1j1 *= NumMatrixTypes;
#include <cstring>

#include <cstring>
    // compute forward scores
#include <cstring>
    for (int i = 0; i <= seq1Length; i++){
#include <cstring>
      unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
#include <cstring>
      for (int j = 0; j <= seq2Length; j++){
#include <cstring>
        unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];
#include <cstring>

#include <cstring>
	if (i > 1 || j > 1){
#include <cstring>
	  if (i > 0 && j > 0){
#include <cstring>
	    forward[0 + ij] = forward[0 + i1j1] + transProb[0][0];
#include <cstring>
	    for (int k = 1; k < NumMatrixTypes; k++)
#include <cstring>
	      LOG_PLUS_EQUALS (forward[0 + ij], forward[k + i1j1] + transProb[k][0]);
#include <cstring>
	    forward[0 + ij] += matchProb[c1][c2];
#include <cstring>
	  }
#include <cstring>
	  if (i > 0){
#include <cstring>
	    for (int k = 0; k < NumInsertStates; k++)
#include <cstring>
	      forward[2*k+1 + ij] = insProb[c1][k] +
#include <cstring>
		LOG_ADD (forward[0 + i1j] + transProb[0][2*k+1],
#include <cstring>
			 forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1]);
#include <cstring>
	  }
#include <cstring>
	  if (j > 0){
#include <cstring>
	    for (int k = 0; k < NumInsertStates; k++)
#include <cstring>
	      forward[2*k+2 + ij] = insProb[c2][k] +
#include <cstring>
		LOG_ADD (forward[0 + ij1] + transProb[0][2*k+2],
#include <cstring>
			 forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2]);
#include <cstring>
	  }
#include <cstring>
	}
#include <cstring>

#include <cstring>
        ij += NumMatrixTypes;
#include <cstring>
        i1j += NumMatrixTypes;
#include <cstring>
        ij1 += NumMatrixTypes;
#include <cstring>
        i1j1 += NumMatrixTypes;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    return forwardPtr;
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeBackwardMatrix()
#include <cstring>
  //
#include <cstring>
  // Computes a set of backward probability matrices for aligning
#include <cstring>
  // seq1 and seq2.
#include <cstring>
  //
#include <cstring>
  // For efficiency reasons, a single-dimensional floating-point
#include <cstring>
  // array is used here, with the following indexing scheme:
#include <cstring>
  //
#include <cstring>
  //    backward[i + NumMatrixTypes * (j * (seq2Length+1) + k)]
#include <cstring>
  //    refers to the probability of starting in state i and
#include <cstring>
  //    aligning from character j+1 to the end of the first
#include <cstring>
  //    sequence and from character k+1 to the end of the second
#include <cstring>
  //    sequence.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  VF *ComputeBackwardMatrix (Sequence *seq1, Sequence *seq2) const {
#include <cstring>

#include <cstring>
    assert (seq1);
#include <cstring>
    assert (seq2);
#include <cstring>

#include <cstring>
    const int seq1Length = seq1->GetLength();
#include <cstring>
    const int seq2Length = seq2->GetLength();
#include <cstring>
    SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
#include <cstring>
    SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
#include <cstring>

#include <cstring>
    // create matrix
#include <cstring>
    VF *backwardPtr = new VF (NumMatrixTypes * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
#include <cstring>
    assert (backwardPtr);
#include <cstring>
    VF &backward = *backwardPtr;
#include <cstring>

#include <cstring>
    // initialization condition
#include <cstring>
    for (int k = 0; k < NumMatrixTypes; k++)
#include <cstring>
      backward[NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1) + k] = initialDistribution[k];
#include <cstring>

#include <cstring>
    // remember offset for each index combination
#include <cstring>
    int ij = (seq1Length+1) * (seq2Length+1) - 1;
#include <cstring>
    int i1j = ij + seq2Length + 1;
#include <cstring>
    int ij1 = ij + 1;
#include <cstring>
    int i1j1 = ij + seq2Length + 2;
#include <cstring>

#include <cstring>
    ij *= NumMatrixTypes;
#include <cstring>
    i1j *= NumMatrixTypes;
#include <cstring>
    ij1 *= NumMatrixTypes;
#include <cstring>
    i1j1 *= NumMatrixTypes;
#include <cstring>

#include <cstring>
    // compute backward scores
#include <cstring>
    for (int i = seq1Length; i >= 0; i--){
#include <cstring>
      unsigned char c1 = (i == seq1Length) ? '~' : (unsigned char) iter1[i+1];
#include <cstring>
      for (int j = seq2Length; j >= 0; j--){
#include <cstring>
        unsigned char c2 = (j == seq2Length) ? '~' : (unsigned char) iter2[j+1];
#include <cstring>

#include <cstring>
        if (i < seq1Length && j < seq2Length){
#include <cstring>
          const float ProbXY = backward[0 + i1j1] + matchProb[c1][c2];
#include <cstring>
          for (int k = 0; k < NumMatrixTypes; k++)
#include <cstring>
            LOG_PLUS_EQUALS (backward[k + ij], ProbXY + transProb[k][0]);
#include <cstring>
        }
#include <cstring>
        if (i < seq1Length){
#include <cstring>
          for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
            LOG_PLUS_EQUALS (backward[0 + ij], backward[2*k+1 + i1j] + insProb[c1][k] + transProb[0][2*k+1]);
#include <cstring>
            LOG_PLUS_EQUALS (backward[2*k+1 + ij], backward[2*k+1 + i1j] + insProb[c1][k] + transProb[2*k+1][2*k+1]);
#include <cstring>
          }
#include <cstring>
        }
#include <cstring>
        if (j < seq2Length){
#include <cstring>
          for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
            LOG_PLUS_EQUALS (backward[0 + ij], backward[2*k+2 + ij1] + insProb[c2][k] + transProb[0][2*k+2]);
#include <cstring>
            LOG_PLUS_EQUALS (backward[2*k+2 + ij], backward[2*k+2 + ij1] + insProb[c2][k] + transProb[2*k+2][2*k+2]);
#include <cstring>
          }
#include <cstring>
        }
#include <cstring>

#include <cstring>
        ij -= NumMatrixTypes;
#include <cstring>
        i1j -= NumMatrixTypes;
#include <cstring>
        ij1 -= NumMatrixTypes;
#include <cstring>
        i1j1 -= NumMatrixTypes;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    return backwardPtr;
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeTotalProbability()
#include <cstring>
  //
#include <cstring>
  // Computes the total probability of an alignment given
#include <cstring>
  // the forward and backward matrices.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  float ComputeTotalProbability (int seq1Length, int seq2Length,
#include <cstring>
                                 const VF &forward, const VF &backward) const {
#include <cstring>

#include <cstring>
    // compute total probability
#include <cstring>
    float totalForwardProb = LOG_ZERO;
#include <cstring>
    float totalBackwardProb = LOG_ZERO;
#include <cstring>
    for (int k = 0; k < NumMatrixTypes; k++){
#include <cstring>
      LOG_PLUS_EQUALS (totalForwardProb,
#include <cstring>
                       forward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] + 
#include <cstring>
		       backward[k + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
#include <cstring>
    }
#include <cstring>

#include <cstring>
    totalBackwardProb = 
#include <cstring>
      forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] +
#include <cstring>
      backward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)];
#include <cstring>

#include <cstring>
    for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
      LOG_PLUS_EQUALS (totalBackwardProb,
#include <cstring>
		       forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] +
#include <cstring>
		       backward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)]);
#include <cstring>
      LOG_PLUS_EQUALS (totalBackwardProb,
#include <cstring>
		       forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] +
#include <cstring>
		       backward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)]);
#include <cstring>
    }
#include <cstring>

#include <cstring>
    //    cerr << totalForwardProb << " " << totalBackwardProb << endl;
#include <cstring>
    
#include <cstring>
    return (totalForwardProb + totalBackwardProb) / 2;
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputePosteriorMatrix()
#include <cstring>
  //
#include <cstring>
  // Computes the posterior probability matrix based on
#include <cstring>
  // the forward and backward matrices.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  VF *ComputePosteriorMatrix (Sequence *seq1, Sequence *seq2,
#include <cstring>
                              const VF &forward, const VF &backward) const {
#include <cstring>

#include <cstring>
    assert (seq1);
#include <cstring>
    assert (seq2);
#include <cstring>

#include <cstring>
    const int seq1Length = seq1->GetLength();
#include <cstring>
    const int seq2Length = seq2->GetLength();
#include <cstring>

#include <cstring>
    float totalProb = ComputeTotalProbability (seq1Length, seq2Length,
#include <cstring>
                                               forward, backward);
#include <cstring>

#include <cstring>
    // compute posterior matrices
#include <cstring>
    VF *posteriorPtr = new VF((seq1Length+1) * (seq2Length+1)); assert (posteriorPtr);
#include <cstring>
    VF &posterior = *posteriorPtr;
#include <cstring>

#include <cstring>
    int ij = 0;
#include <cstring>
    VF::iterator ptr = posterior.begin();
#include <cstring>

#include <cstring>
    for (int i = 0; i <= seq1Length; i++){
#include <cstring>
      for (int j = 0; j <= seq2Length; j++){
#include <cstring>
        *(ptr++) = EXP (min (LOG_ONE, forward[ij] + backward[ij] - totalProb));
#include <cstring>
        ij += NumMatrixTypes;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    posterior[0] = 0;
#include <cstring>

#include <cstring>
    return posteriorPtr;
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /*
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeExpectedCounts()
#include <cstring>
  //
#include <cstring>
  // Computes the expected counts for the various transitions.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  VVF *ComputeExpectedCounts () const {
#include <cstring>

#include <cstring>
    assert (seq1);
#include <cstring>
    assert (seq2);
#include <cstring>

#include <cstring>
    const int seq1Length = seq1->GetLength();
#include <cstring>
    const int seq2Length = seq2->GetLength();
#include <cstring>
    SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
#include <cstring>
    SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
#include <cstring>

#include <cstring>
    // compute total probability
#include <cstring>
    float totalProb = ComputeTotalProbability (seq1Length, seq2Length,
#include <cstring>
                                               forward, backward);
#include <cstring>

#include <cstring>
    // initialize expected counts
#include <cstring>
    VVF *countsPtr = new VVF(NumMatrixTypes + 1, VF(NumMatrixTypes, LOG_ZERO)); assert (countsPtr);
#include <cstring>
    VVF &counts = *countsPtr;
#include <cstring>

#include <cstring>
    // remember offset for each index combination
#include <cstring>
    int ij = 0;
#include <cstring>
    int i1j = -seq2Length - 1;
#include <cstring>
    int ij1 = -1;
#include <cstring>
    int i1j1 = -seq2Length - 2;
#include <cstring>

#include <cstring>
    ij *= NumMatrixTypes;
#include <cstring>
    i1j *= NumMatrixTypes;
#include <cstring>
    ij1 *= NumMatrixTypes;
#include <cstring>
    i1j1 *= NumMatrixTypes;
#include <cstring>

#include <cstring>
    // compute expected counts
#include <cstring>
    for (int i = 0; i <= seq1Length; i++){
#include <cstring>
      unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
#include <cstring>
      for (int j = 0; j <= seq2Length; j++){
#include <cstring>
        unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];
#include <cstring>

#include <cstring>
        if (i > 0 && j > 0){
#include <cstring>
          for (int k = 0; k < NumMatrixTypes; k++)
#include <cstring>
            LOG_PLUS_EQUALS (counts[k][0],
#include <cstring>
                             forward[k + i1j1] + transProb[k][0] +
#include <cstring>
                             matchProb[c1][c2] + backward[0 + ij]);
#include <cstring>
        }
#include <cstring>
        if (i > 0){
#include <cstring>
          for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
            LOG_PLUS_EQUALS (counts[0][2*k+1],
#include <cstring>
                             forward[0 + i1j] + transProb[0][2*k+1] +
#include <cstring>
                             insProb[c1][k] + backward[2*k+1 + ij]);
#include <cstring>
            LOG_PLUS_EQUALS (counts[2*k+1][2*k+1],
#include <cstring>
                             forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1] +
#include <cstring>
                             insProb[c1][k] + backward[2*k+1 + ij]);
#include <cstring>
          }
#include <cstring>
        }
#include <cstring>
        if (j > 0){
#include <cstring>
          for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
            LOG_PLUS_EQUALS (counts[0][2*k+2],
#include <cstring>
                             forward[0 + ij1] + transProb[0][2*k+2] +
#include <cstring>
                             insProb[c2][k] + backward[2*k+2 + ij]);
#include <cstring>
            LOG_PLUS_EQUALS (counts[2*k+2][2*k+2],
#include <cstring>
                             forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2] +
#include <cstring>
                             insProb[c2][k] + backward[2*k+2 + ij]);
#include <cstring>
          }
#include <cstring>
        }
#include <cstring>

#include <cstring>
        ij += NumMatrixTypes;
#include <cstring>
        i1j += NumMatrixTypes;
#include <cstring>
        ij1 += NumMatrixTypes;
#include <cstring>
        i1j1 += NumMatrixTypes;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // scale all expected counts appropriately
#include <cstring>
    for (int i = 0; i < NumMatrixTypes; i++)
#include <cstring>
      for (int j = 0; j < NumMatrixTypes; j++)
#include <cstring>
        counts[i][j] -= totalProb;
#include <cstring>

#include <cstring>
  }
#include <cstring>
  */
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeNewParameters()
#include <cstring>
  //
#include <cstring>
  // Computes a new parameter set based on the expected counts
#include <cstring>
  // given.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  void ComputeNewParameters (Sequence *seq1, Sequence *seq2,
#include <cstring>
			     const VF &forward, const VF &backward,
#include <cstring>
                             VF &initDistribMat, VF &gapOpen,
#include <cstring>
                             VF &gapExtend, VVF &emitPairs, VF &emitSingle, bool enableTrainEmissions) const {
#include <cstring>
    
#include <cstring>
    assert (seq1);
#include <cstring>
    assert (seq2);
#include <cstring>

#include <cstring>
    const int seq1Length = seq1->GetLength();
#include <cstring>
    const int seq2Length = seq2->GetLength();
#include <cstring>
    SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
#include <cstring>
    SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
#include <cstring>

#include <cstring>
    // compute total probability
#include <cstring>
    float totalProb = ComputeTotalProbability (seq1Length, seq2Length,
#include <cstring>
                                               forward, backward);
#include <cstring>
    
#include <cstring>
    // initialize expected counts
#include <cstring>
    VVF transCounts (NumMatrixTypes, VF (NumMatrixTypes, LOG_ZERO));
#include <cstring>
    VF initCounts (NumMatrixTypes, LOG_ZERO);
#include <cstring>
    VVF pairCounts (256, VF (256, LOG_ZERO));
#include <cstring>
    VF singleCounts (256, LOG_ZERO);
#include <cstring>
    
#include <cstring>
    // remember offset for each index combination
#include <cstring>
    int ij = 0;
#include <cstring>
    int i1j = -seq2Length - 1;
#include <cstring>
    int ij1 = -1;
#include <cstring>
    int i1j1 = -seq2Length - 2;
#include <cstring>

#include <cstring>
    ij *= NumMatrixTypes;
#include <cstring>
    i1j *= NumMatrixTypes;
#include <cstring>
    ij1 *= NumMatrixTypes;
#include <cstring>
    i1j1 *= NumMatrixTypes;
#include <cstring>

#include <cstring>
    // compute initial distribution posteriors
#include <cstring>
    initCounts[0] = LOG_ADD (forward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)] +
#include <cstring>
			     backward[0 + NumMatrixTypes * (1 * (seq2Length+1) + 1)],
#include <cstring>
			     forward[0 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] + 
#include <cstring>
			     backward[0 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
#include <cstring>
    for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
      initCounts[2*k+1] = LOG_ADD (forward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)] +
#include <cstring>
				   backward[2*k+1 + NumMatrixTypes * (1 * (seq2Length+1) + 0)],
#include <cstring>
				   forward[2*k+1 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] + 
#include <cstring>
				   backward[2*k+1 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
#include <cstring>
      initCounts[2*k+2] = LOG_ADD (forward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)] +
#include <cstring>
				   backward[2*k+2 + NumMatrixTypes * (0 * (seq2Length+1) + 1)],
#include <cstring>
				   forward[2*k+2 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)] + 
#include <cstring>
				   backward[2*k+2 + NumMatrixTypes * ((seq1Length+1) * (seq2Length+1) - 1)]);
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // compute expected counts
#include <cstring>
    for (int i = 0; i <= seq1Length; i++){
#include <cstring>
      unsigned char c1 = (i == 0) ? '~' : (unsigned char) toupper(iter1[i]);
#include <cstring>
      for (int j = 0; j <= seq2Length; j++){
#include <cstring>
        unsigned char c2 = (j == 0) ? '~' : (unsigned char) toupper(iter2[j]);
#include <cstring>

#include <cstring>
	if (i > 0 && j > 0){
#include <cstring>
	  if (enableTrainEmissions && i == 1 && j == 1){
#include <cstring>
	    LOG_PLUS_EQUALS (pairCounts[c1][c2],
#include <cstring>
			     initialDistribution[0] + matchProb[c1][c2] + backward[0 + ij]);
#include <cstring>
	    LOG_PLUS_EQUALS (pairCounts[c2][c1],
#include <cstring>
			     initialDistribution[0] + matchProb[c2][c1] + backward[0 + ij]);
#include <cstring>
	  }
#include <cstring>

#include <cstring>
	  for (int k = 0; k < NumMatrixTypes; k++){
#include <cstring>
	    LOG_PLUS_EQUALS (transCounts[k][0],
#include <cstring>
			     forward[k + i1j1] + transProb[k][0] +
#include <cstring>
			     matchProb[c1][c2] + backward[0 + ij]);
#include <cstring>
	    if (enableTrainEmissions && i != 1 || j != 1){
#include <cstring>
	      LOG_PLUS_EQUALS (pairCounts[c1][c2],
#include <cstring>
			       forward[k + i1j1] + transProb[k][0] +
#include <cstring>
			       matchProb[c1][c2] + backward[0 + ij]);
#include <cstring>
	      LOG_PLUS_EQUALS (pairCounts[c2][c1],
#include <cstring>
			       forward[k + i1j1] + transProb[k][0] +
#include <cstring>
			       matchProb[c2][c1] + backward[0 + ij]);
#include <cstring>
	    }
#include <cstring>
	  }
#include <cstring>
	}
#include <cstring>
	if (i > 0){
#include <cstring>
	  for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
	    LOG_PLUS_EQUALS (transCounts[0][2*k+1],
#include <cstring>
			     forward[0 + i1j] + transProb[0][2*k+1] +
#include <cstring>
			     insProb[c1][k] + backward[2*k+1 + ij]);
#include <cstring>
	    LOG_PLUS_EQUALS (transCounts[2*k+1][2*k+1],
#include <cstring>
			     forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1] +
#include <cstring>
			     insProb[c1][k] + backward[2*k+1 + ij]);
#include <cstring>
	    if (enableTrainEmissions){
#include <cstring>
	      if (i == 1 && j == 0){
#include <cstring>
		LOG_PLUS_EQUALS (singleCounts[c1],
#include <cstring>
				 initialDistribution[2*k+1] + insProb[c1][k] + backward[2*k+1 + ij]);
#include <cstring>
	      }
#include <cstring>
	      else {
#include <cstring>
		LOG_PLUS_EQUALS (singleCounts[c1],
#include <cstring>
				 forward[0 + i1j] + transProb[0][2*k+1] +
#include <cstring>
				 insProb[c1][k] + backward[2*k+1 + ij]);
#include <cstring>
		LOG_PLUS_EQUALS (singleCounts[c1],
#include <cstring>
				 forward[2*k+1 + i1j] + transProb[2*k+1][2*k+1] +
#include <cstring>
				 insProb[c1][k] + backward[2*k+1 + ij]);
#include <cstring>
	      }
#include <cstring>
	    }
#include <cstring>
	  }
#include <cstring>
	}
#include <cstring>
	if (j > 0){
#include <cstring>
	  for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
	    LOG_PLUS_EQUALS (transCounts[0][2*k+2],
#include <cstring>
			     forward[0 + ij1] + transProb[0][2*k+2] +
#include <cstring>
			     insProb[c2][k] + backward[2*k+2 + ij]);
#include <cstring>
	    LOG_PLUS_EQUALS (transCounts[2*k+2][2*k+2],
#include <cstring>
			     forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2] +
#include <cstring>
			     insProb[c2][k] + backward[2*k+2 + ij]);
#include <cstring>
	    if (enableTrainEmissions){
#include <cstring>
	      if (i == 0 && j == 1){
#include <cstring>
		LOG_PLUS_EQUALS (singleCounts[c2],
#include <cstring>
				 initialDistribution[2*k+2] + insProb[c2][k] + backward[2*k+2 + ij]);
#include <cstring>
	      }
#include <cstring>
	      else {
#include <cstring>
		LOG_PLUS_EQUALS (singleCounts[c2],
#include <cstring>
				 forward[0 + ij1] + transProb[0][2*k+2] +
#include <cstring>
				 insProb[c2][k] + backward[2*k+2 + ij]);
#include <cstring>
		LOG_PLUS_EQUALS (singleCounts[c2],
#include <cstring>
				 forward[2*k+2 + ij1] + transProb[2*k+2][2*k+2] +
#include <cstring>
				 insProb[c2][k] + backward[2*k+2 + ij]);
#include <cstring>
	      }
#include <cstring>
	    }
#include <cstring>
	  }
#include <cstring>
	}
#include <cstring>
      
#include <cstring>
        ij += NumMatrixTypes;
#include <cstring>
        i1j += NumMatrixTypes;
#include <cstring>
        ij1 += NumMatrixTypes;
#include <cstring>
        i1j1 += NumMatrixTypes;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // scale all expected counts appropriately
#include <cstring>
    for (int i = 0; i < NumMatrixTypes; i++){
#include <cstring>
      initCounts[i] -= totalProb;
#include <cstring>
      for (int j = 0; j < NumMatrixTypes; j++)
#include <cstring>
        transCounts[i][j] -= totalProb;
#include <cstring>
    }
#include <cstring>
    if (enableTrainEmissions){
#include <cstring>
      for (int i = 0; i < 256; i++){
#include <cstring>
	for (int j = 0; j < 256; j++)
#include <cstring>
	  pairCounts[i][j] -= totalProb;
#include <cstring>
	singleCounts[i] -= totalProb;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // compute new initial distribution
#include <cstring>
    float totalInitDistribCounts = 0;
#include <cstring>
    for (int i = 0; i < NumMatrixTypes; i++)
#include <cstring>
      totalInitDistribCounts += exp (initCounts[i]); // should be 2
#include <cstring>
    initDistribMat[0] = min (1.0f, max (0.0f, (float) exp (initCounts[0]) / totalInitDistribCounts));
#include <cstring>
    for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
      float val = (exp (initCounts[2*k+1]) + exp (initCounts[2*k+2])) / 2;
#include <cstring>
      initDistribMat[2*k+1] = initDistribMat[2*k+2] = min (1.0f, max (0.0f, val / totalInitDistribCounts));
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // compute total counts for match state
#include <cstring>
    float inMatchStateCounts = 0;
#include <cstring>
    for (int i = 0; i < NumMatrixTypes; i++)
#include <cstring>
      inMatchStateCounts += exp (transCounts[0][i]);
#include <cstring>
    for (int i = 0; i < NumInsertStates; i++){
#include <cstring>

#include <cstring>
      // compute total counts for gap state
#include <cstring>
      float inGapStateCounts =
#include <cstring>
        exp (transCounts[2*i+1][0]) +
#include <cstring>
        exp (transCounts[2*i+1][2*i+1]) +
#include <cstring>
        exp (transCounts[2*i+2][0]) +
#include <cstring>
        exp (transCounts[2*i+2][2*i+2]);
#include <cstring>

#include <cstring>
      gapOpen[2*i] = gapOpen[2*i+1] =
#include <cstring>
        (exp (transCounts[0][2*i+1]) +
#include <cstring>
         exp (transCounts[0][2*i+2])) /
#include <cstring>
        (2 * inMatchStateCounts);
#include <cstring>

#include <cstring>
      gapExtend[2*i] = gapExtend[2*i+1] =
#include <cstring>
        (exp (transCounts[2*i+1][2*i+1]) +
#include <cstring>
         exp (transCounts[2*i+2][2*i+2])) /
#include <cstring>
        inGapStateCounts;
#include <cstring>
    }
#include <cstring>

#include <cstring>
    if (enableTrainEmissions){
#include <cstring>
      float totalPairCounts = 0;
#include <cstring>
      float totalSingleCounts = 0;
#include <cstring>
      for (int i = 0; i < 256; i++){
#include <cstring>
	for (int j = 0; j <= i; j++)
#include <cstring>
	  totalPairCounts += exp (pairCounts[j][i]);
#include <cstring>
	totalSingleCounts += exp (singleCounts[i]);
#include <cstring>
      }
#include <cstring>
      
#include <cstring>
      for (int i = 0; i < 256; i++) if (!islower ((char) i)){
#include <cstring>
	int li = (int)((unsigned char) tolower ((char) i));
#include <cstring>
	for (int j = 0; j <= i; j++) if (!islower ((char) j)){
#include <cstring>
	  int lj = (int)((unsigned char) tolower ((char) j));
#include <cstring>
	  emitPairs[i][j] = emitPairs[i][lj] = emitPairs[li][j] = emitPairs[li][lj] = 
#include <cstring>
	    emitPairs[j][i] = emitPairs[j][li] = emitPairs[lj][i] = emitPairs[lj][li] = exp(pairCounts[j][i]) / totalPairCounts;
#include <cstring>
	}
#include <cstring>
	emitSingle[i] = emitSingle[li] = exp(singleCounts[i]) / totalSingleCounts;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>
  }
#include <cstring>
    
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeAlignment()
#include <cstring>
  //
#include <cstring>
  // Computes an alignment based on the given posterior matrix.
#include <cstring>
  // This is done by finding the maximum summing path (or
#include <cstring>
  // maximum weight trace) through the posterior matrix.  The
#include <cstring>
  // final alignment is returned as a pair consisting of:
#include <cstring>
  //    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
#include <cstring>
  //        denote insertions in one of the two sequences and
#include <cstring>
  //        B's denote that both sequences are present (i.e.
#include <cstring>
  //        matches).
#include <cstring>
  //    (2) a float indicating the sum achieved
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  pair<SafeVector<char> *, float> ComputeAlignment (int seq1Length, int seq2Length,
#include <cstring>
                                                    const VF &posterior) const {
#include <cstring>

#include <cstring>
    float *twoRows = new float[(seq2Length+1)*2]; assert (twoRows);
#include <cstring>
    float *oldRow = twoRows;
#include <cstring>
    float *newRow = twoRows + seq2Length + 1;
#include <cstring>

#include <cstring>
    char *tracebackMatrix = new char[(seq1Length+1)*(seq2Length+1)]; assert (tracebackMatrix);
#include <cstring>
    char *tracebackPtr = tracebackMatrix;
#include <cstring>

#include <cstring>
    VF::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;
#include <cstring>

#include <cstring>
    // initialization
#include <cstring>
    for (int i = 0; i <= seq2Length; i++){
#include <cstring>
      oldRow[i] = 0;
#include <cstring>
      *(tracebackPtr++) = 'L';
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // fill in matrix
#include <cstring>
    for (int i = 1; i <= seq1Length; i++){
#include <cstring>

#include <cstring>
      // initialize left column
#include <cstring>
      newRow[0] = 0;
#include <cstring>
      posteriorPtr++;
#include <cstring>
      *(tracebackPtr++) = 'U';
#include <cstring>

#include <cstring>
      // fill in rest of row
#include <cstring>
      for (int j = 1; j <= seq2Length; j++){
#include <cstring>
        ChooseBestOfThree (*(posteriorPtr++) + oldRow[j-1], newRow[j-1], oldRow[j],
#include <cstring>
                           'D', 'L', 'U', &newRow[j], tracebackPtr++);
#include <cstring>
      }
#include <cstring>

#include <cstring>
      // swap rows
#include <cstring>
      float *temp = oldRow;
#include <cstring>
      oldRow = newRow;
#include <cstring>
      newRow = temp;
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // store best score
#include <cstring>
    float total = oldRow[seq2Length];
#include <cstring>
    delete [] twoRows;
#include <cstring>

#include <cstring>
    // compute traceback
#include <cstring>
    SafeVector<char> *alignment = new SafeVector<char>; assert (alignment);
#include <cstring>
    int r = seq1Length, c = seq2Length;
#include <cstring>
    while (r != 0 || c != 0){
#include <cstring>
      char ch = tracebackMatrix[r*(seq2Length+1) + c];
#include <cstring>
      switch (ch){
#include <cstring>
      case 'L': c--; alignment->push_back ('Y'); break;
#include <cstring>
      case 'U': r--; alignment->push_back ('X'); break;
#include <cstring>
      case 'D': c--; r--; alignment->push_back ('B'); break;
#include <cstring>
      default: assert (false);
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    delete [] tracebackMatrix;
#include <cstring>

#include <cstring>
    reverse (alignment->begin(), alignment->end());
#include <cstring>

#include <cstring>
    return make_pair(alignment, total);
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeAlignmentWithGapPenalties()
#include <cstring>
  //
#include <cstring>
  // Similar to ComputeAlignment() except with gap penalties.
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  pair<SafeVector<char> *, float> ComputeAlignmentWithGapPenalties (MultiSequence *align1,
#include <cstring>
                                                                    MultiSequence *align2,
#include <cstring>
                                                                    const VF &posterior, int numSeqs1,
#include <cstring>
                                                                    int numSeqs2,
#include <cstring>
                                                                    float gapOpenPenalty,
#include <cstring>
                                                                    float gapContinuePenalty) const {
#include <cstring>
    int seq1Length = align1->GetSequence(0)->GetLength();
#include <cstring>
    int seq2Length = align2->GetSequence(0)->GetLength();
#include <cstring>
    SafeVector<SafeVector<char>::iterator > dataPtrs1 (align1->GetNumSequences());
#include <cstring>
    SafeVector<SafeVector<char>::iterator > dataPtrs2 (align2->GetNumSequences());
#include <cstring>

#include <cstring>
    // grab character data
#include <cstring>
    for (int i = 0; i < align1->GetNumSequences(); i++)
#include <cstring>
      dataPtrs1[i] = align1->GetSequence(i)->GetDataPtr();
#include <cstring>
    for (int i = 0; i < align2->GetNumSequences(); i++)
#include <cstring>
      dataPtrs2[i] = align2->GetSequence(i)->GetDataPtr();
#include <cstring>

#include <cstring>
    // the number of active sequences at any given column is defined to be the
#include <cstring>
    // number of non-gap characters in that column; the number of gap opens at
#include <cstring>
    // any given column is defined to be the number of gap characters in that
#include <cstring>
    // column where the previous character in the respective sequence was not
#include <cstring>
    // a gap
#include <cstring>
    SafeVector<int> numActive1 (seq1Length+1), numGapOpens1 (seq1Length+1);
#include <cstring>
    SafeVector<int> numActive2 (seq2Length+1), numGapOpens2 (seq2Length+1);
#include <cstring>

#include <cstring>
    // compute number of active sequences and gap opens for each group
#include <cstring>
    for (int i = 0; i < align1->GetNumSequences(); i++){
#include <cstring>
      SafeVector<char>::iterator dataPtr = align1->GetSequence(i)->GetDataPtr();
#include <cstring>
      numActive1[0] = numGapOpens1[0] = 0;
#include <cstring>
      for (int j = 1; j <= seq1Length; j++){
#include <cstring>
        if (dataPtr[j] != '-'){
#include <cstring>
          numActive1[j]++;
#include <cstring>
          numGapOpens1[j] += (j != 1 && dataPtr[j-1] != '-');
#include <cstring>
        }
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>
    for (int i = 0; i < align2->GetNumSequences(); i++){
#include <cstring>
      SafeVector<char>::iterator dataPtr = align2->GetSequence(i)->GetDataPtr();
#include <cstring>
      numActive2[0] = numGapOpens2[0] = 0;
#include <cstring>
      for (int j = 1; j <= seq2Length; j++){
#include <cstring>
        if (dataPtr[j] != '-'){
#include <cstring>
          numActive2[j]++;
#include <cstring>
          numGapOpens2[j] += (j != 1 && dataPtr[j-1] != '-');
#include <cstring>
        }
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    VVF openingPenalty1 (numSeqs1+1, VF (numSeqs2+1));
#include <cstring>
    VF continuingPenalty1 (numSeqs1+1);
#include <cstring>
    VVF openingPenalty2 (numSeqs1+1, VF (numSeqs2+1));
#include <cstring>
    VF continuingPenalty2 (numSeqs2+1);
#include <cstring>

#include <cstring>
    // precompute penalties
#include <cstring>
    for (int i = 0; i <= numSeqs1; i++)
#include <cstring>
      for (int j = 0; j <= numSeqs2; j++)
#include <cstring>
        openingPenalty1[i][j] = i * (gapOpenPenalty * j + gapContinuePenalty * (numSeqs2 - j));
#include <cstring>
    for (int i = 0; i <= numSeqs1; i++)
#include <cstring>
      continuingPenalty1[i] = i * gapContinuePenalty * numSeqs2;
#include <cstring>
    for (int i = 0; i <= numSeqs2; i++)
#include <cstring>
      for (int j = 0; j <= numSeqs1; j++)
#include <cstring>
        openingPenalty2[i][j] = i * (gapOpenPenalty * j + gapContinuePenalty * (numSeqs1 - j));
#include <cstring>
    for (int i = 0; i <= numSeqs2; i++)
#include <cstring>
      continuingPenalty2[i] = i * gapContinuePenalty * numSeqs1;
#include <cstring>

#include <cstring>
    float *twoRows = new float[6*(seq2Length+1)]; assert (twoRows);
#include <cstring>
    float *oldRowMatch = twoRows;
#include <cstring>
    float *newRowMatch = twoRows + (seq2Length+1);
#include <cstring>
    float *oldRowInsertX = twoRows + 2*(seq2Length+1);
#include <cstring>
    float *newRowInsertX = twoRows + 3*(seq2Length+1);
#include <cstring>
    float *oldRowInsertY = twoRows + 4*(seq2Length+1);
#include <cstring>
    float *newRowInsertY = twoRows + 5*(seq2Length+1);
#include <cstring>

#include <cstring>
    char *tracebackMatrix = new char[3*(seq1Length+1)*(seq2Length+1)]; assert (tracebackMatrix);
#include <cstring>
    char *tracebackPtr = tracebackMatrix;
#include <cstring>

#include <cstring>
    VF::const_iterator posteriorPtr = posterior.begin() + seq2Length + 1;
#include <cstring>

#include <cstring>
    // initialization
#include <cstring>
    for (int i = 0; i <= seq2Length; i++){
#include <cstring>
      oldRowMatch[i] = oldRowInsertX[i] = (i == 0) ? 0 : LOG_ZERO;
#include <cstring>
      oldRowInsertY[i] = (i == 0) ? 0 : oldRowInsertY[i-1] + continuingPenalty2[numActive2[i]];
#include <cstring>
      *(tracebackPtr) = *(tracebackPtr+1) = *(tracebackPtr+2) = 'Y';
#include <cstring>
      tracebackPtr += 3;
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // fill in matrix
#include <cstring>
    for (int i = 1; i <= seq1Length; i++){
#include <cstring>

#include <cstring>
      // initialize left column
#include <cstring>
      newRowMatch[0] = newRowInsertY[0] = LOG_ZERO;
#include <cstring>
      newRowInsertX[0] = oldRowInsertX[0] + continuingPenalty1[numActive1[i]];
#include <cstring>
      posteriorPtr++;
#include <cstring>
      *(tracebackPtr) = *(tracebackPtr+1) = *(tracebackPtr+2) = 'X';
#include <cstring>
      tracebackPtr += 3;
#include <cstring>

#include <cstring>
      // fill in rest of row
#include <cstring>
      for (int j = 1; j <= seq2Length; j++){
#include <cstring>

#include <cstring>
        // going to MATCH state
#include <cstring>
        ChooseBestOfThree (oldRowMatch[j-1],
#include <cstring>
                           oldRowInsertX[j-1],
#include <cstring>
                           oldRowInsertY[j-1],
#include <cstring>
                           'M', 'X', 'Y', &newRowMatch[j], tracebackPtr++);
#include <cstring>
        newRowMatch[j] += *(posteriorPtr++);
#include <cstring>

#include <cstring>
        // going to INSERT X state
#include <cstring>
        ChooseBestOfThree (oldRowMatch[j] + openingPenalty1[numActive1[i]][numGapOpens2[j]],
#include <cstring>
                           oldRowInsertX[j] + continuingPenalty1[numActive1[i]],
#include <cstring>
                           oldRowInsertY[j] + openingPenalty1[numActive1[i]][numGapOpens2[j]],
#include <cstring>
                           'M', 'X', 'Y', &newRowInsertX[j], tracebackPtr++);
#include <cstring>

#include <cstring>
        // going to INSERT Y state
#include <cstring>
        ChooseBestOfThree (newRowMatch[j-1] + openingPenalty2[numActive2[j]][numGapOpens1[i]],
#include <cstring>
                           newRowInsertX[j-1] + openingPenalty2[numActive2[j]][numGapOpens1[i]],
#include <cstring>
                           newRowInsertY[j-1] + continuingPenalty2[numActive2[j]],
#include <cstring>
                           'M', 'X', 'Y', &newRowInsertY[j], tracebackPtr++);
#include <cstring>
      }
#include <cstring>

#include <cstring>
      // swap rows
#include <cstring>
      float *temp;
#include <cstring>
      temp = oldRowMatch; oldRowMatch = newRowMatch; newRowMatch = temp;
#include <cstring>
      temp = oldRowInsertX; oldRowInsertX = newRowInsertX; newRowInsertX = temp;
#include <cstring>
      temp = oldRowInsertY; oldRowInsertY = newRowInsertY; newRowInsertY = temp;
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // store best score
#include <cstring>
    float total;
#include <cstring>
    char matrix;
#include <cstring>
    ChooseBestOfThree (oldRowMatch[seq2Length], oldRowInsertX[seq2Length], oldRowInsertY[seq2Length],
#include <cstring>
                       'M', 'X', 'Y', &total, &matrix);
#include <cstring>

#include <cstring>
    delete [] twoRows;
#include <cstring>

#include <cstring>
    // compute traceback
#include <cstring>
    SafeVector<char> *alignment = new SafeVector<char>; assert (alignment);
#include <cstring>
    int r = seq1Length, c = seq2Length;
#include <cstring>
    while (r != 0 || c != 0){
#include <cstring>

#include <cstring>
      int offset = (matrix == 'M') ? 0 : (matrix == 'X') ? 1 : 2;
#include <cstring>
      char ch = tracebackMatrix[(r*(seq2Length+1) + c) * 3 + offset];
#include <cstring>
      switch (matrix){
#include <cstring>
      case 'Y': c--; alignment->push_back ('Y'); break;
#include <cstring>
      case 'X': r--; alignment->push_back ('X'); break;
#include <cstring>
      case 'M': c--; r--; alignment->push_back ('B'); break;
#include <cstring>
      default: assert (false);
#include <cstring>
      }
#include <cstring>
      matrix = ch;
#include <cstring>
    }
#include <cstring>

#include <cstring>
    delete [] tracebackMatrix;
#include <cstring>

#include <cstring>
    reverse (alignment->begin(), alignment->end());
#include <cstring>

#include <cstring>
    return make_pair(alignment, 1.0f);
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::ComputeViterbiAlignment()
#include <cstring>
  //
#include <cstring>
  // Computes the highest probability pairwise alignment using the
#include <cstring>
  // probabilistic model.  The final alignment is returned as a
#include <cstring>
  //  pair consisting of:
#include <cstring>
  //    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
#include <cstring>
  //        denote insertions in one of the two sequences and
#include <cstring>
  //        B's denote that both sequences are present (i.e.
#include <cstring>
  //        matches).
#include <cstring>
  //    (2) a float containing the log probability of the best
#include <cstring>
  //        alignment (not used)
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  pair<SafeVector<char> *, float> ComputeViterbiAlignment (Sequence *seq1, Sequence *seq2) const {
#include <cstring>
    
#include <cstring>
    assert (seq1);
#include <cstring>
    assert (seq2);
#include <cstring>
    
#include <cstring>
    const int seq1Length = seq1->GetLength();
#include <cstring>
    const int seq2Length = seq2->GetLength();
#include <cstring>
    
#include <cstring>
    // retrieve the points to the beginning of each sequence
#include <cstring>
    SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
#include <cstring>
    SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
#include <cstring>
    
#include <cstring>
    // create viterbi matrix
#include <cstring>
    VF *viterbiPtr = new VF (NumMatrixTypes * (seq1Length+1) * (seq2Length+1), LOG_ZERO);
#include <cstring>
    assert (viterbiPtr);
#include <cstring>
    VF &viterbi = *viterbiPtr;
#include <cstring>

#include <cstring>
    // create traceback matrix
#include <cstring>
    VI *tracebackPtr = new VI (NumMatrixTypes * (seq1Length+1) * (seq2Length+1), -1);
#include <cstring>
    assert (tracebackPtr);
#include <cstring>
    VI &traceback = *tracebackPtr;
#include <cstring>

#include <cstring>
    // initialization condition
#include <cstring>
    for (int k = 0; k < NumMatrixTypes; k++)
#include <cstring>
      viterbi[k] = initialDistribution[k];
#include <cstring>

#include <cstring>
    // remember offset for each index combination
#include <cstring>
    int ij = 0;
#include <cstring>
    int i1j = -seq2Length - 1;
#include <cstring>
    int ij1 = -1;
#include <cstring>
    int i1j1 = -seq2Length - 2;
#include <cstring>

#include <cstring>
    ij *= NumMatrixTypes;
#include <cstring>
    i1j *= NumMatrixTypes;
#include <cstring>
    ij1 *= NumMatrixTypes;
#include <cstring>
    i1j1 *= NumMatrixTypes;
#include <cstring>

#include <cstring>
    // compute viterbi scores
#include <cstring>
    for (int i = 0; i <= seq1Length; i++){
#include <cstring>
      unsigned char c1 = (i == 0) ? '~' : (unsigned char) iter1[i];
#include <cstring>
      for (int j = 0; j <= seq2Length; j++){
#include <cstring>
        unsigned char c2 = (j == 0) ? '~' : (unsigned char) iter2[j];
#include <cstring>

#include <cstring>
        if (i > 0 && j > 0){
#include <cstring>
          for (int k = 0; k < NumMatrixTypes; k++){
#include <cstring>
	    float newVal = viterbi[k + i1j1] + transProb[k][0] + matchProb[c1][c2];
#include <cstring>
	    if (viterbi[0 + ij] < newVal){
#include <cstring>
	      viterbi[0 + ij] = newVal;
#include <cstring>
	      traceback[0 + ij] = k;
#include <cstring>
	    }
#include <cstring>
	  }
#include <cstring>
        }
#include <cstring>
        if (i > 0){
#include <cstring>
          for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
	    float valFromMatch = insProb[c1][k] + viterbi[0 + i1j] + transProb[0][2*k+1];
#include <cstring>
	    float valFromIns = insProb[c1][k] + viterbi[2*k+1 + i1j] + transProb[2*k+1][2*k+1];
#include <cstring>
	    if (valFromMatch >= valFromIns){
#include <cstring>
	      viterbi[2*k+1 + ij] = valFromMatch;
#include <cstring>
	      traceback[2*k+1 + ij] = 0;
#include <cstring>
	    }
#include <cstring>
	    else {
#include <cstring>
	      viterbi[2*k+1 + ij] = valFromIns;
#include <cstring>
	      traceback[2*k+1 + ij] = 2*k+1;
#include <cstring>
	    }
#include <cstring>
	  }
#include <cstring>
	}
#include <cstring>
        if (j > 0){
#include <cstring>
          for (int k = 0; k < NumInsertStates; k++){
#include <cstring>
	    float valFromMatch = insProb[c2][k] + viterbi[0 + ij1] + transProb[0][2*k+2];
#include <cstring>
	    float valFromIns = insProb[c2][k] + viterbi[2*k+2 + ij1] + transProb[2*k+2][2*k+2];
#include <cstring>
	    if (valFromMatch >= valFromIns){
#include <cstring>
	      viterbi[2*k+2 + ij] = valFromMatch;
#include <cstring>
	      traceback[2*k+2 + ij] = 0;
#include <cstring>
	    }
#include <cstring>
	    else {
#include <cstring>
	      viterbi[2*k+2 + ij] = valFromIns;
#include <cstring>
	      traceback[2*k+2 + ij] = 2*k+2;
#include <cstring>
	    }
#include <cstring>
	  }
#include <cstring>
        }
#include <cstring>

#include <cstring>
        ij += NumMatrixTypes;
#include <cstring>
        i1j += NumMatrixTypes;
#include <cstring>
        ij1 += NumMatrixTypes;
#include <cstring>
        i1j1 += NumMatrixTypes;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>

#include <cstring>
    // figure out best terminating cell
#include <cstring>
    float bestProb = LOG_ZERO;
#include <cstring>
    int state = -1;
#include <cstring>
    for (int k = 0; k < NumMatrixTypes; k++){
#include <cstring>
      float thisProb = viterbi[k + NumMatrixTypes * ((seq1Length+1)*(seq2Length+1) - 1)] + initialDistribution[k];
#include <cstring>
      if (bestProb < thisProb){
#include <cstring>
	bestProb = thisProb;
#include <cstring>
	state = k;
#include <cstring>
      }
#include <cstring>
    }
#include <cstring>
    assert (state != -1);
#include <cstring>

#include <cstring>
    delete viterbiPtr;
#include <cstring>

#include <cstring>
    // compute traceback
#include <cstring>
    SafeVector<char> *alignment = new SafeVector<char>; assert (alignment);
#include <cstring>
    int r = seq1Length, c = seq2Length;
#include <cstring>
    while (r != 0 || c != 0){
#include <cstring>
      int newState = traceback[state + NumMatrixTypes * (r * (seq2Length+1) + c)];
#include <cstring>
      
#include <cstring>
      if (state == 0){ c--; r--; alignment->push_back ('B'); }
#include <cstring>
      else if (state % 2 == 1){ r--; alignment->push_back ('X'); }
#include <cstring>
      else { c--; alignment->push_back ('Y'); }
#include <cstring>
      
#include <cstring>
      state = newState;
#include <cstring>
    }
#include <cstring>

#include <cstring>
    delete tracebackPtr;
#include <cstring>

#include <cstring>
    reverse (alignment->begin(), alignment->end());
#include <cstring>
    
#include <cstring>
    return make_pair(alignment, bestProb);
#include <cstring>
  }
#include <cstring>

#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>
  // ProbabilisticModel::BuildPosterior()
#include <cstring>
  //
#include <cstring>
  // Builds a posterior probability matrix needed to align a pair
#include <cstring>
  // of alignments.  Mathematically, the returned matrix M is
#include <cstring>
  // defined as follows:
#include <cstring>
  //    M[i,j] =     sum          sum      f(s,t,i,j)
#include <cstring>
  //             s in align1  t in align2
#include <cstring>
  // where
#include <cstring>
  //                  [  P(s[i'] <--> t[j'])
#include <cstring>
  //                  [       if s[i'] is a letter in the ith column of align1 and
#include <cstring>
  //                  [          t[j'] it a letter in the jth column of align2
#include <cstring>
  //    f(s,t,i,j) =  [
#include <cstring>
  //                  [  0    otherwise
#include <cstring>
  //
#include <cstring>
  /////////////////////////////////////////////////////////////////
#include <cstring>

#include <cstring>
  VF *BuildPosterior (MultiSequence *align1, MultiSequence *align2,
#include <cstring>
                      const SafeVector<SafeVector<SparseMatrix *> > &sparseMatrices,
#include <cstring>
		      float cutoff = 0.0f) const {
#include <cstring>
    const int seq1Length = align1->GetSequence(0)->GetLength();
#include <cstring>
    const int seq2Length = align2->GetSequence(0)->GetLength();
#include <cstring>

#include <cstring>
    VF *posteriorPtr = new VF((seq1Length+1) * (seq2Length+1), 0); assert (posteriorPtr);
#include <cstring>
    VF &posterior = *posteriorPtr;
#include <cstring>
    VF::iterator postPtr = posterior.begin();
#include <cstring>

#include <cstring>
    // for each s in align1
#include <cstring>
    for (int i = 0; i < align1->GetNumSequences(); i++){
#include <cstring>
      int first = align1->GetSequence(i)->GetLabel();
#include <cstring>
      SafeVector<int> *mapping1 = align1->GetSequence(i)->GetMapping();
#include <cstring>

#include <cstring>
      // for each t in align2
#include <cstring>
      for (int j = 0; j < align2->GetNumSequences(); j++){
#include <cstring>
        int second = align2->GetSequence(j)->GetLabel();
#include <cstring>
        SafeVector<int> *mapping2 = align2->GetSequence(j)->GetMapping();
#include <cstring>

#include <cstring>
	if (first < second){
#include <cstring>

#include <cstring>
	  // get the associated sparse matrix
#include <cstring>
	  SparseMatrix *matrix = sparseMatrices[first][second];
#include <cstring>
	  
#include <cstring>
	  for (int ii = 1; ii <= matrix->GetSeq1Length(); ii++){
#include <cstring>
	    SafeVector<PIF>::iterator row = matrix->GetRowPtr(ii);
#include <cstring>
	    int base = (*mapping1)[ii] * (seq2Length+1);
#include <cstring>
	    int rowSize = matrix->GetRowSize(ii);
#include <cstring>
	    
#include <cstring>
	    // add in all relevant values
#include <cstring>
	    for (int jj = 0; jj < rowSize; jj++)
#include <cstring>
	      posterior[base + (*mapping2)[row[jj].first]] += row[jj].second;
#include <cstring>
	    
#include <cstring>
	    // subtract cutoff 
#include <cstring>
	    for (int jj = 0; jj < matrix->GetSeq2Length(); jj++)
#include <cstring>
	      posterior[base + (*mapping2)[jj]] -= cutoff;
#include <cstring>
	  }
#include <cstring>

#include <cstring>
	} else {
#include <cstring>

#include <cstring>
	  // get the associated sparse matrix
#include <cstring>
	  SparseMatrix *matrix = sparseMatrices[second][first];
#include <cstring>
	  
#include <cstring>
	  for (int jj = 1; jj <= matrix->GetSeq1Length(); jj++){
#include <cstring>
	    SafeVector<PIF>::iterator row = matrix->GetRowPtr(jj);
#include <cstring>
	    int base = (*mapping2)[jj];
#include <cstring>
	    int rowSize = matrix->GetRowSize(jj);
#include <cstring>
	    
#include <cstring>
	    // add in all relevant values
#include <cstring>
	    for (int ii = 0; ii < rowSize; ii++)
#include <cstring>
	      posterior[base + (*mapping1)[row[ii].first] * (seq2Length + 1)] += row[ii].second;
#include <cstring>
	    
#include <cstring>
	    // subtract cutoff 
#include <cstring>
	    for (int ii = 0; ii < matrix->GetSeq2Length(); ii++)
#include <cstring>
	      posterior[base + (*mapping1)[ii] * (seq2Length + 1)] -= cutoff;
#include <cstring>
	  }
#include <cstring>

#include <cstring>
	}
#include <cstring>
	
#include <cstring>

#include <cstring>
        delete mapping2;
#include <cstring>
      }
#include <cstring>

#include <cstring>
      delete mapping1;
#include <cstring>
    }
#include <cstring>

#include <cstring>
    return posteriorPtr;
#include <cstring>
  }
#include <cstring>
};
#include <cstring>

#include <cstring>
#endif
#include <cstring>
