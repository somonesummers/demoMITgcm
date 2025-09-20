Modified version of Iceplume by Tom Cowton.
The only code modification is the multiplication of the melt rate by 1/sin_beta in iceplume_cal.F
sin_beta is sin(beta) and beta is the angle between the horizontal and the glacier face. If the glacier face is vertical sin_beta=1.

As specified in Bonneau et al (2023), to implement a constant slope, the user has to multiplied gravity (g) by (1/sin_beta)^2 and
the entrainment (E_0) by 1/sin_beta. 

If someone wish to use this modification, they only have to replace the file in the pkg folder of the iceplume_package by these ones. 
