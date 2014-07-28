#!/usr/bin/perl -w

# This perl script is used to invoke stilts (if it's installed and the "stilts"
# command properly calls it) to generate a plot of retained mass fraction for
# the different orbits generated by the SALTSA_example.cpp script.

$cmd = "stilts plot2d ofmt=eps out='SALTSA_example_orbits.eps' legend=false " .
	"grid=false ylog=false xpix=400 ypix=400 ylo=-0.8 ";

for($i=0;$i<=9;$i++) {
    $cmd = $cmd . "ifmt".$i."=ascii in".$i."=SALTSA_example_orbit_".$i.".dat " .
	"cmd".$i."='addcol log\\(M/M0\\) log10(m_ret)' ydata".$i."=log\\(M/M0\\) " .
	"cmd".$i."='addcol t/P t/8.77683' xdata".$i."=t/P " .
        "shape".$i."=filled_circle " .
	"size".$i."=0 line".$i."=DotToDot "
}

for($i=0;$i<=2;$i++) {
    $cmd = $cmd . "ifmt1".$i."=ascii in1".$i."=SALTSA_example_orbit_1".$i.".dat " .
	"cmd1".$i."='addcol log\\(M/M0\\) log10(m_ret)' ydata1".$i."=log\\(M/M0\\) " .
	"cmd1".$i."='addcol t/P t/8.77683' xdata1".$i."=t/P " .
        "shape1".$i."=filled_circle " .
	"size1".$i."=0 line1".$i."=DotToDot "
}

$cmd = $cmd . "color0=black color1=blue color2=cyan color3=green color4=yellow " .
              "color5=red color6=magenta color7=blue color8=cyan color9=green " .
              "color10=yellow color11=red color12=black ";

$cmd = $cmd . "name0='c=1.00' name1='c=.99' name2='c=.95' name3='c=.9' " .
    "name4='c=.8' name5='c=.7' name6='c=.6' name7='c=.5' name8='c=.4' " .
    "name9='c=.3' name10='c=.2' name11='c=.1' name12='c=.05'";

system($cmd);
