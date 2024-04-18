The program for performance and stability calculation of various rotor craft - RotaryWingSim
The version number is 2.x to distinguish it from the pure helicopter program 1.0. 
The program is a beta release which means it is far from complete or thoroughly tested, 
rather anyone using it is a beta tester and I would appreciate your contributions to enhance the program. 
Throughout the four years of working on this program I have made numerous attempts to find people supporting my effort, 
so far without success. One major reason for releasing the program thus is one more attempt to find groups or individuals 
who will help to improve the program. You can support this goal by spreading word about this program even if you do not 
intend to use it or take part in its further development yourself. Please notify anyone who might be interested in this 
program by talking to students, universities or whoever you deem likely to be inclined to participate in this effort.

Program Features
Currently the program offers two basic calculation modes. A helicopter mode where the rotational speed of the rotor is 
constant and an autogyro mode where either rotor rpm for a given collective pitch or collective pitch for a given rotor 
rpm are calculated. The different calculation modes currently offer the features listed below ( prescribed input values 
are to the left of the -> , calculated output is to the right).

a) Helicopter Mode
- Performance Calculation: Rate of Climb, Power to Hover, Flight Envelope
- Trim Calculation:
- Prescribed Rotational Speed -> Cyclic, Collective Pitch

b) Autogyro Mode
- Trim Calculation 1:
...Prescribed Rotational Speed -> Cyclic and Collective Pitch
- Trim Calculation 2:
...Prescribed Collective Pitch -> Rotational Speed, Cyclic Pitch
-
c) Helicopter- and Autogyro Mode
- Calculation of aircraft trim values for straight and level as well as
...climbing/descending flight.
- Calculation of Engine power to trim.
- Linear geometric blade twist can be included
- Rotor state calculation method (a): formulae of naca 716 for rectangular
...blades are used. The blades are considered rigid, no cyclic twist
...considered
- Rotor state calculation method (b): formulae of naca 600 for rectangular
...blades are used. Cyclic twist as exhibited e.g. by the Cierva C-30 is included
- Especially for helicopters the horizontal stabilizer is sometimes connected
...to the cyclic via a linkage that changes the angle of attack in conjunction
...with cyclic input. Such so called gearing curves can be entered as
...polynomials.
- Rotor downwash on the horizontal stabilizer and the fuselage is included
...using the values from Heyson and Katzoff
- Fuselage aerodynamic coefficients are included using polynomials
- Investigation of winged or wingless configurations.
- Wing aerodynamic coefficients are included using two dimensional,
...Reynolds dependent profile data in xfoil format.
- Stability: stability derivatives, solution of stability quartic

In both modes the aircraft flight state in trimmed flight is calculated first, then the stability derivatives are 
calculated via small perturbations applied to each of the trim variables in turn. The stability derivatives in 
autogyro mode are currently not deemed correct since rotor speed does not change during the perturbation steps used 
to calculate them.


The archive contains the code, sample input files and a brief user guide.

Version 2p1 is the one which allows calculation of helicopter and autogyros. The simulink model in this file does no longer work. 
It is included to offer a basis from which to work if you want one.
You will need Octave, a free Matlab (TM) clone, or Matlab(R) to run the program. 
Please read the "Getting Started" section of the user guide
