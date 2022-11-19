#include "../include/amjuel_hyd_recombination.hxx"

/// Coefficients to calculate the effective reaction rate <σv>
/// Reaction 2.1.8, Amjuel page 141 (section H.4)
/// E-index varies fastest, so coefficient is [T][n]
static constexpr const BoutReal rate_coefs[9][9] = {
    {-28.58858570847, 0.02068671746773, -0.007868331504755, 0.003843362133859,
     -0.0007411492158905, 9.273687892997e-05, -7.063529824805e-06, 3.026539277057e-07,
     -5.373940838104e-09},
    {-0.7676413320499, 0.0127800603259, -0.01870326896978, 0.00382855504889,
     -0.0003627770385335, 4.401007253801e-07, 1.932701779173e-06, -1.176872895577e-07,
     2.215851843121e-09},
    {0.002823851790251, -0.001907812518731, 0.01121251125171, -0.003711328186517,
     0.0006617485083301, -6.860774445002e-05, 4.508046989099e-06, -1.723423509284e-07,
     2.805361431741e-09},
    {-0.01062884273731, -0.01010719783828, 0.004208412930611, -0.00100574441054,
     0.0001013652422369, -2.044691594727e-06, -4.431181498017e-07, 3.457903389784e-08,
     -7.374639775683e-10},
    {0.001582701550903, 0.002794099401979, -0.002024796037098, 0.0006250304936976,
     -9.224891301052e-05, 7.546853961575e-06, -3.682709551169e-07, 1.035928615391e-08,
     -1.325312585168e-10},
    {-0.0001938012790522, 0.0002148453735781, 3.393285358049e-05, -3.746423753955e-05,
     7.509176112468e-06, -8.688365258514e-07, 7.144767938783e-08, -3.367897014044e-09,
     6.250111099227e-11},
    {6.041794354114e-06, -0.0001421502819671, 6.14387907608e-05, -1.232549226121e-05,
     1.394562183496e-06, -6.434833988001e-08, -2.746804724917e-09, 3.564291012995e-10,
     -8.55170819761e-12},
    {1.742316850715e-06, 1.595051038326e-05, -7.858419208668e-06, 1.774935420144e-06,
     -2.187584251561e-07, 1.327090702659e-08, -1.386720240985e-10, -1.946206688519e-11,
     5.745422385081e-13},
    {-1.384927774988e-07, -5.664673433879e-07, 2.886857762387e-07, -6.591743182569e-08,
     8.008790343319e-09, -4.805837071646e-10, 6.459706573699e-12, 5.510729582791e-13,
     -1.680871303639e-14}};

/// Coefficients to calculate the radiation energy loss
/// Reaction 2.1.8, Amjuel page 284 (section H.10)
/// E-index varies fastest, so coefficient is [T][n]
static constexpr const BoutReal radiation_coefs[9][9] = {
    {-25.92450349909, 0.01222097271874, 4.278499401907e-05, 0.001943967743593,
     -0.0007123474602102, 0.0001303523395892, -1.186560752561e-05, 5.334455630031e-07,
     -9.349857887253e-09},
    {-0.7290670236493, -0.01540323930666, -0.00340609377919, 0.001532243431817,
     -0.0004658423772784, 5.972448753445e-05, -4.070843294052e-06, 1.378709880644e-07,
     -1.818079729166e-09},
    {0.02363925869096, 0.01164453346305, -0.005845209334594, 0.002854145868307,
     -0.0005077485291132, 4.211106637742e-05, -1.251436618314e-06, -1.626555745259e-08,
     1.073458810743e-09},
    {0.003645333930947, -0.001005820792983, 0.0006956352274249, -0.0009305056373739,
     0.0002584896294384, -3.294643898894e-05, 2.112924018518e-06, -6.544682842175e-08,
     7.8102930757e-10},
    {0.001594184648757, -1.582238007548e-05, 0.0004073695619272, -9.379169243859e-05,
     1.490890502214e-06, 2.245292872209e-06, -3.150901014513e-07, 1.631965635818e-08,
     -2.984093025695e-10},
    {-0.001216668033378, -0.0003503070140126, 0.0001043500296633, 9.536162767321e-06,
     -6.908681884097e-06, 8.232019008169e-07, -2.905331051259e-08, -3.169038517749e-10,
     2.442765766167e-11},
    {0.0002376115895241, 0.0001172709777146, -6.695182045674e-05, 1.18818400621e-05,
     -4.381514364966e-07, -6.936267173079e-08, 6.592249255001e-09, -1.778887958831e-10,
     1.160762106747e-12},
    {-1.930977636766e-05, -1.318401491304e-05, 8.848025453481e-06, -2.07237071139e-06,
     2.055919993599e-07, -7.489632654212e-09, -7.073797030749e-11, 1.047087505147e-11,
     -1.87744627135e-13},
    {5.599257775146e-07, 4.977823319311e-07, -3.615013823092e-07, 9.466989306497e-08,
     -1.146485227699e-08, 6.772338917155e-10, -1.776496344763e-11, 7.199195061382e-14,
     3.929300283002e-15}};

void AmjuelHydRecombination::calculate_rates(
  Options& electron, Options& atom, Options& ion, 
  Field3D &reaction_rate, Field3D &momentum_exchange,
  Field3D &energy_exchange, Field3D &energy_loss) {
  electron_reaction(electron, ion, atom, rate_coefs, radiation_coefs,
                    13.6, // Potential energy loss [eV] heats electrons
                    reaction_rate, momentum_exchange, energy_exchange, energy_loss
  );
}
