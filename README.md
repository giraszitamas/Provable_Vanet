Provably secure authenticated anonymous batch messaging for VANETs
Andrea Huszti, Norbert Olah

Vehicles require the ability to communicate with each other to improve transportation efficiency
and safety for both vehicles and pedestrians. To address this, the Intelligent Transport System recommends
applying Vehicle Ad-Hoc Networks (VANETs), which use cryptographic protocols to enable authorized
vehicles to report road conditions anonymously. Our previous solution ([1]), which uses identity-based
cryptography, has been improved by reducing the number of operations required by the On-Board Unit
(OBU) during incident reporting, resulting in increased efficiency. By leveraging bilinear pairing properties,
the master secret key need not be stored, and our recommendation enables batch verification of messages
and revocation of sender anonymity. We introduce a new adversarial model and a definition for a secure
anonymous authenticated message broadcast scheme and show that our scheme is secure if Computational
Diffie-Hellman problem is computationally infeasible, and the bilinear map is one-way, where values of the
map are taken from the Anonymous User List. Additionally, the implementation of the scheme demonstrates
the good efficiency results
