/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file simulation_cec.hpp
  \brief Simulation-based CEC

  EPFL CS-472 2021 Final Project Option 2
*/

#pragma once

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operations.hpp>

#include "../utils/node_map.hpp"
#include "miter.hpp"
#include "simulation.hpp"


namespace mockturtle
{

/* Statistics to be reported */
struct simulation_cec_stats
{
  /*! \brief Split variable (simulation size). */
  uint32_t split_var{ 0 };

  /*! \brief Number of simulation rounds. */
  uint32_t rounds{ 0 };
};

namespace detail
{

template<class Ntk>
class simulation_cec_impl
{
public:
  using pattern_t = unordered_node_map<kitty::dynamic_truth_table, Ntk>;
  using node = typename Ntk::node;
  using signal = typename Ntk::signal;

public:
  explicit simulation_cec_impl( Ntk& ntk, simulation_cec_stats& st )
      : _ntk( ntk ),
        _st( st )
  {
  }

  bool run()
  {
    /* TODO: write your implementation here */

    /*compute split_var and rounds */
    _st.split_var = compute_splitting_var(_ntk);
    _st.rounds = compute_rounds(_ntk,_st.split_var);

    /*print split_var and rounds*/
    std::cout << "split_var : " << _st.split_var << std::endl;
    std::cout << "num of rounds : " << _st.rounds << std::endl;

    pattern_t patterns(_ntk);


    /*initialize  patterns and simulate first round */ 
    init_patterns( _ntk, _st.split_var,patterns);
    default_simulator<kitty::dynamic_truth_table> simulator(_st.split_var);
    simulate_nodes( _ntk, patterns,simulator );

    /*first equivalence check, returning only False if the equivalence check is False */
    if (!check_eq( _ntk,patterns )){
      return(false);
    }
   
    /*update and simulate the different rounds */
    for (uint32_t round = 1; round <= _st.rounds - 1; ++round){

       pattern_clear(patterns);
       update_pattern(round,patterns);
       simulate_nodes( _ntk, patterns,simulator );

      /* quit the loop only if a check is False, if it is True we need to iterate through all the simulation rounds*/
        if (!check_eq( _ntk,patterns )){
          return(false);
          }
     }
    
    return true;
    
  }

private:
  /* you can add additional methods here */

  /* compute split_var based on the mathematical expression  */
  uint32_t compute_splitting_var(Ntk& _N){
    uint32_t n,v,k;
    n=_N.num_pis();   /* number of inputs */
    v=_N._storage->nodes.size();  /*number of nodes */
    k= 3+ ( (log((1<<29)/v-32))/log(2) ) ; /*upper limit of m */
    if (n<=6) {
      return(n);
    } 
    else {
      return(std::min(k,n));   /* split_var should satisfy both conditions <n and <=k , hence the function min(k,n)  */
    }
  }

  /* compute number of rounds of the simulation  */
  uint32_t compute_rounds(Ntk& _N, uint32_t split_var){
    uint32_t n=_N.num_pis();
    return(1 << (n-split_var));
  }

  void init_patterns( Ntk& _N, uint32_t split_var, pattern_t& patterns){

     /*pattern_t patterns(_N);*/
     
     _N.foreach_pi( [&]( auto const& n, auto p ){ 
      kitty::dynamic_truth_table tt (split_var);
      if (p < split_var) {
      kitty:: create_nth_var(tt , p);
      }
      patterns[n]=tt;
    } );

  }
  
/* equivalence check function */
  bool check_eq (Ntk& _N ,pattern_t& patterns){
    bool equiv_check ;
    equiv_check=true;
    _N.foreach_po( [&]( auto const& f) {
    if ( _N.is_complemented( f ) )  /* the case of complemented function*/
    {
      if ( !is_const0(~patterns[f]) ) {  /*checking bit by bit if the outputs are 0 (True) */
        equiv_check = false;
      }
    }
    else
    {
      if ( !is_const0(patterns[f]) ) {
      equiv_check = false; 
      }
    }
  } );
  return equiv_check;
}

/*the function to update the pattern*/
  void update_pattern( uint32_t& round ,pattern_t& patterns ){
    uint32_t r = round;  

      _ntk.foreach_pi( [&]( auto const& j, auto k ) 
      {
        if (k >= _st.split_var ){
          if (r % 2 == 1) {
             if ( is_const0(patterns[j]) ) patterns[j] = ~patterns[j];
          }
      
        else {
           if ( !is_const0(patterns[j]) ) patterns[j] = ~patterns[j];
        }
        r /= 2;
  
        }
      } );
  }

    /* function to clear the patterns to avoid alloc errors */
     void pattern_clear( pattern_t& patterns ){
    _ntk.foreach_gate( [&]( auto const& n )
    {
       patterns.erase(n);
    } );
  }

//};

  
private:
  Ntk& _ntk;
  simulation_cec_stats& _st;
  /* you can add other attributes here */
};

} // namespace detail

/* Entry point for users to call */

/*! \brief Simulation-based CEC.
 *
 * This function implements a simulation-based combinational equivalence checker.
 * The implementation creates a miter network and run several rounds of simulation
 * to verify the functional equivalence. For memory and speed reasons this approach
 * is limited up to 40 input networks. It returns an optional which is `nullopt`,
 * if the network has more than 40 inputs.
 */
template<class Ntk>
std::optional<bool> simulation_cec( Ntk const& ntk1, Ntk const& ntk2, simulation_cec_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_num_pis_v<Ntk>, "Ntk does not implement the num_pis method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_pi_v<Ntk>, "Ntk does not implement the foreach_pi method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );

  simulation_cec_stats st;

  bool result = false;

  if ( ntk1.num_pis() > 40 )
    return std::nullopt;

  auto ntk_miter = miter<Ntk>( ntk1, ntk2 );

  if ( ntk_miter.has_value() )
  {
    detail::simulation_cec_impl p( *ntk_miter, st );
    result = p.run();
  }

  if ( pst )
    *pst = st;

  return result;
}

} // namespace mockturtle
