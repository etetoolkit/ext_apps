/******************************COPYRIGHT NOTICE*******************************/
/*  (c) Centro de Regulacio Genomica                                                        */
/*  and                                                                                     */
/*  Cedric Notredame                                                                        */
/*  22 Oct 2013 - 15:51.                                                                    */
/*All rights reserved.                                                                      */
/*This file is part of T-COFFEE.                                                            */
/*                                                                                          */
/*    T-COFFEE is free software; you can redistribute it and/or modify                      */
/*    it under the terms of the GNU General Public License as published by                  */
/*    the Free Software Foundation; either version 2 of the License, or                     */
/*    (at your option) any later version.                                                   */
/*                                                                                          */
/*    T-COFFEE is distributed in the hope that it will be useful,                           */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of                        */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                         */
/*    GNU General Public License for more details.                                          */
/*                                                                                          */
/*    You should have received a copy of the GNU General Public License                     */
/*    along with Foobar; if not, write to the Free Software                                 */
/*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA             */
/*...............................................                                           */
/*  If you need some more information                                                       */
/*  cedric.notredame@europe.com                                                             */
/*...............................................                                           */
/******************************COPYRIGHT NOTICE*******************************/
typedef struct 
{
  int N_COMPONENT;
  double *double_logB_alpha;
  double *exponant_list;
  double **ALPHA;
  double *DM_Q;
  double *alpha_tot;
  int n_aa;
  int tot_n;
}
Mixture;

double int_logB (int *i, int n);
double float_logB (float *i, int n);
double double_logB (double *x, int n);
double *** make_lup_table ( Mixture *D);
double  double_logB2(int j, double *n,Mixture *D);
double compute_exponant ( double *n, int j, Mixture *D);

double *compute_matrix_p ( double *n,int Nseq);
double* compute_dirichlet_p ( double *n,int Nseq);
void precompute_log_B ( double *table,Mixture *D);
double compute_X (double *n,int i,Mixture *D);
Mixture *read_dirichlet ( char *name);
int dirichlet_code( char aa);

double lgamma2 ( double x);
double lgamma_r(double x, int *signgamp);
