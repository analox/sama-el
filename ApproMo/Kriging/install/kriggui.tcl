#!/usr/bin/wish
#krigify.tcl
# Anthony Padula, Rice University

global version
set version {Version 1.1, Last Updated 7/12/2000}
global InitialValue         ;# initial & final values of editable entries
global FinalValue          ;# in Tcl/Tk dialogs
# 0 = seed, 1 = stream, 3 = savefile, 4 = loadfile
# 5,6 convert
set FinalValue(0) -1
set FinalValue(1) 42
set FinalValue(3) ""
set FinalValue(4) ""
set FinalValue(5) ""
set FinalValue(6) ""
global _ind
global oldval
global WTitle
global Windows
#exec xset b on
set soundflag 1 
# 0 == make sound, 1 == silent
#
# plotting == 0 when active, 1 when waiting for plot to be calculated
set plotting 0
set barwid 10
set bgcol white
set status {Status: Set parameters and push Plot}
set fname {Picture: }
#set seedstring {Your Seed:            Your Stream:}
set p 2
set alpha 2
# need a slider for
set n 50
set tau2 1
set beta0 0
set theta 1
set sigma2 100
set lower 0
set upper 10

set WTitle "Krigifier"
set Windows ""

#font create deffont -family Courier -size 18 -weight normal
font create deffont -size 18

# make a noise if appropriate
proc noise {} { 
    global soundflag
    global myturn
#puts stdout "soundflag = $soundflag"
    if {$soundflag == 0} { 
	bell
    }
}

# display current status in appropriate window
proc displaystatus {} {
    global status
    .cs itemconfigure stattag -text $status
}

proc displayfname {} {
    global fname
    .cs itemconfigure ftag -text $fname
}

proc displayseed {} {
    global seedstring
    .cs itemconfigure stag -text $seedstring
}

# to process a move
proc replot {} {
    global plotting
    global status
    global p
    global n
    global tau2
    global beta0
    global alpha
    global theta
    global sigma2
    global lower
    global upper
    global FinalValue
    if {$plotting == 0} {
	#send move to client
	puts stdout "$p $n $tau2 $beta0 $alpha $theta $sigma2 $lower $upper\
		$FinalValue(0) $FinalValue(1)"
	flush stdout
	#update and display current status
	set status {Status: Plotting...}
	set plotting 1
    } else { }
    displaystatus
}

proc changesound {} {
  global soundflag;
  if {$soundflag == 0} {
    set soundflag 1
    .bf.sound configure -text "Sound"
  } else {
    set soundflag 0
    .bf.sound configure -text "Silent"
  }
}

proc changep {} {
  global p
  if {$p == 1} {
    set p 2
  } else {
    set p 1
  }
  .bf.pbut configure -text "p = $p"
}

proc changealpha {} {
  global alpha
  if {$alpha == 1} {
    set alpha 2
  } else {
    set alpha 1
  }
  .bf.alphabut configure -text "alpha = $alpha"
}

proc changelower {a} {
  global lower
  global upper
  if {[expr "$lower" - "$upper"] >= 0} {
    noise
    set lower [expr "$upper" - 1]
  }
}

proc changeupper {a} {
  global lower
  global upper
  if {[expr "$lower" - "$upper"] >= 0} {
    noise
    set upper [expr "$lower" + 1]
  }
}


#  Here we have several nice procedures borrowed from Bill Bynum's
#  submit script http://www.cs.wm.edu/~bynum/classtools/install_submit

######################################################################
#
#  Window handling procs
#     necessary for a graceful exit from X
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  add_window
#     adds window 'w' to the Window list, if it's not already there
proc add_window { w } {
   global Windows
   set i [lsearch -exact $Windows $w]
   if { $i < 0 } {   ;# only add 'w' if it's not there
      lappend Windows $w
   }
}  ;# add_window


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  destroy_window
#     destoys window 'w'
#     and if it was the only window left, kills the withdrawn .
proc destroy_window { w } {
   global Windows
   set i [lsearch -exact $Windows $w]
   if { $i >= 0 } {
      set Windows [lreplace $Windows $i $i]
   }
   destroy $w
#   if { [llength $Windows] == 0} {
 #     destroy .
  # }
}  ;# destroy_window


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  make_entry_frame     
#     creates an editable entry frame for a Tk display having the form     
#                          
#                              +---------------------------------+
#                              |             +=================+ |
# frame, named "frame_name"    | label_text  |   init_text     | |
#                              |             +=================+ |
#                              +---------------------------------+
#                              |<----------->|<----------------->|
#                                 lt_width       it_width 
#
# in the window "wn"
#
# adds the two symbols 
#    "InitialValue($valindex)" and "FinalValue($valindex)"
# to the two global arrays used by the caller to  
# compare the initial value of the editable entry text and its
# value after editing
#
proc make_entry_frame {wn frame_name valindex label_text lt_width init_text it_width} {
   global InitialValue
   global FinalValue 
   if {$wn == "." } {
      set wn ""
   }
   global $wn.$frame_name
   frame $wn.$frame_name
   label $wn.$frame_name\_l -text "$label_text" -width $lt_width -anchor w
   entry $wn.$frame_name\_e -width $it_width -relief sunken\
      -textvariable FinalValue($valindex)
   set InitialValue($valindex) "$init_text"
   $wn.$frame_name\_e delete 0 end 
   $wn.$frame_name\_e insert end $InitialValue($valindex)
   pack $wn.$frame_name\_l $wn.$frame_name\_e -side left \
      -in $wn.$frame_name
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  ok_win 
#     displays s and waits for click of OK button
#
#  Anthony Padula modified this to include an entry box and cancel button
proc ok_win { s valindex } {
   global WTitle
   global FinalValue
   global _ind
   global oldval

   set _ind $valindex
   set oldval $FinalValue($valindex)
   toplevel .td  
   add_window .td
   set len [string length $s]
   make_entry_frame .td en $valindex "$s" $len "" 30  
   #label .td.l -text "$s"
   button .td.ok -text OK -command { destroy_window .td}
   button .td.cancel -text Cancel -command { set FinalValue($_ind) $oldval; destroy_window .td} 
   wm geometry .td +400+400
   wm title .td "$WTitle -- OK?"
#   pack .td.l -pady 2m -padx 2m
   pack .td.en -side top -fill x -expand 1
   pack .td.ok -fill x -expand 1 -side right
   pack .td.cancel -fill x -expand 1 -side left
   set win [grab current]
   grab set .td
   tkwait window .td
   #if { $win != "" } { grab set $win }
}  ;# ok_win

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  help_win 
#     displays s and waits for click of OK button
#
proc help_win { s } {
   global WTitle
   toplevel .td
   add_window .td
   label .td.l -text "$s" -justify left
   button .td.ok -text OK -command { destroy_window .td}
   wm geometry .td +400+400
   wm title .td "$WTitle -- Help"
   pack .td.l -pady 2m -padx 2m
   pack .td.ok -fill x -expand 1
   #set win [grab current]
   #grab set .td
   tkwait window .td
   #if { $win != "" } { grab set $win }
}  ;# help_win

proc saveparams {} {
    global p
    global n
    global tau2
    global beta0
    global alpha
    global theta
    global sigma2
    global lower
    global upper
    global FinalValue
    global WTitle 
    set WTitle {Krigifier Save}
    ok_win {Save Parameters to file: } 3 
    if [catch [list open $FinalValue(3) "w"] outf] {
	return
    }
    puts $outf "$p\n$n\n$tau2\n$beta0\n$alpha\n$theta\n$sigma2\n$lower\n$upper"
    puts $outf "$FinalValue(0)\n$FinalValue(1)"
    #close $FinalValue(3)
    close $outf
}

proc loadparams {} {
    global p
    global n
    global tau2
    global beta0
    global alpha
    global theta
    global sigma2
    global lower
    global upper
    global FinalValue 
    global WTitle
    set WTitle {Krigifier Load}
    ok_win {Load Parameters from file: } 4 
    if [catch [list open $FinalValue(4) "r+"] infile] {
	return
    }
    gets $infile p
    gets $infile n
    gets $infile tau2
    gets $infile beta0
    gets $infile alpha
    gets $infile theta
    gets $infile sigma2
    gets $infile lower
    gets $infile upper
    gets $infile FinalValue(0)   
    gets $infile FinalValue(1)
    close $infile
}

proc getigcom {infile} {
   # global var
    gets $infile templine
    while {[expr [scan $templine "#%c" temp] > 0]} {
	gets $infile templine
    }
    return [lindex $templine 0]
}
   
#converts data files 
proc convert {flag} {
# flag == 0   gui ->  krigifier
# flag == 1   krigifier -> gui

#use local names
    global FinalValue 
    global WTitle
    set WTitle {Krigifier Convert}
    set FinalValue(5) ""
    ok_win {Convert from file: } 5 
    if [catch [list open $FinalValue(5) "r+"] infile] {
	return
    }
    set FinalValue(6) ""
    ok_win {Convert to file: } 6 
    if [catch [list open $FinalValue(6) "w"] outfile] {
	return
    }

   if {$flag == 0} {

       gets $infile p
       puts $outfile $p
       gets $infile n
       puts $outfile $n
       gets $infile tau2
       gets $infile beta0
       puts $outfile $beta0
# output sufficient zeros
       if {$p == 2} {
	   puts $outfile {0 0}
	   puts $outfile {0 0}
	   puts $outfile {0 0}
	   puts $outfile {0 0}
       } elseif {$p == 1} {
	   puts $outfile {0}
	   puts $outfile {0}
	   puts $outfile {0}
       }
       gets $infile alpha
       puts $outfile $alpha
       gets $infile theta
       puts $outfile $theta
       gets $infile sigma2
       puts $outfile $sigma2
       gets $infile lower
       if {$p == 2} {
	   puts $outfile "$lower $lower"
       } elseif {$p == 1} {
	   puts $outfile $lower
       }
       gets $infile upper
        if {$p == 2} {
	   puts $outfile "$upper $upper"
       } elseif {$p == 1} {
	   puts $outfile $upper
       }
   } else {
       set p [getigcom $infile]
       puts stdout $p
       puts $outfile [lindex $p 0]
       set n [getigcom $infile]
       puts stdout $n
       puts $outfile [lindex $n 0]
       set beta0 [getigcom $infile]
       puts stdout $beta0
       if {$p == 2} {
	   gets $infile crap
	   puts stdout $crap
	   while {[expr [scan $crap "#%c" temp] > 0]} {
	       gets $infile crap
	   }
	   if {[llength $crap] < 2} {
	       gets $infile $crap
	   }
	   gets $infile crap
	   puts stdout $crap
	   while {[expr [scan $crap "#%c" temp] > 0]} {
	       gets $infile crap
	   }
	   puts stdout $crap
	   set i [expr 4 - [llength $crap] ]
	   puts stdout $i
	   while {$i > 0} { 
	       gets $infile crap
	       set i [expr $i - [llength $crap] ]
	   }
	   gets $infile crap
	   puts stdout $crap
	   while {[expr [scan $crap "#%c" temp] > 0]} {
	       gets $infile crap
	   }
	   if {[llength $crap] < 2} {
	       gets $infile $crap
	   } 


       } elseif {$p == 1} {
	   gets $infile crap
	   gets $infile crap
	   gets $infile crap
       }
       puts $outfile 0
       puts $outfile [lindex $beta0 0]
       set alpha [getigcom $infile]
       puts $outfile [lindex $alpha 0]
       set theta [getigcom $infile]
       puts $outfile [lindex $theta 0]
       set sigma2 [getigcom $infile]
       puts $outfile [lindex $sigma2 0]
       set lower [getigcom $infile]
       puts $outfile [lindex $lower 0]
       set upper [getigcom $infile]
       puts $outfile [lindex $upper 0]
       puts $outfile -1
       puts $outfile 42
	
   }
   close $infile 
   close $outfile

}

proc converthelp {} {
    global WTitle
    global helptext
    set WTitle {Krigifier Data Conversion}

set helptext "\tThe gui\ interface\ and\ the\ regular\ krigifier\ class\ utilize
different\ formats\ in\ their\ data\ files.\ \ This\ is\ because\ the\ krigifier\ by
default\ uses\ a\ fixed\ trend\ specified\ by\ four\ parameters.\ \ The\ gui\ by\ default
generates\ a\ random\ trend\ based\ on\ the\ parameters\ tau^2\ and\ beta0.\ \ Further\ the
krigifier\ is\ set\ up\ to\ handle\ some\ basic\ comments\ in\ the\ data\ file, 
while\ the\ gui\ is\ not.\ \ Thus,\ for\ your\ convenience,\ we\ have\ provided\ the
means\ for\ converting\ between\ one\ data\ style\ and\ the\ other.

\tIn\ the\ \"File\"\ menu,\ the\ entry\ \"Gui\ format\ ->\ krigifier\ format\" 
converts\ files\ produced\ by\ the\ gui's\ \"Save\"\ function\ into\ files\ which
can\ be\ read\ by\ the\ krigifier.\ \ In\ order\ to\ reproduce\ a\ function\ you\ have
created\ using\ the\ gui,\ you\ should\ use\ this\ option,\ and\ enter\ the\ names
of\ the\ appropriate\ files\ when\ prompted.\ \ In\ your\ code,\ you\ might\ then\ type:

randfunc\ MyRandomFunction(seed,\ stream);
MyRandomFunction.setvalues(\"krig.dat\");
MyRandomFunction.generateTrend(tau2,\ beta0);

where\ \"seed\",\ \"stream\",\ \"tau2\",\ \"beta0\"\ are\ the\ appropriate\ values
that\ were\ set\ on\ the\ gui,\ and\ \"krig.dat\"\ is\ the\ file\ name\ you\ entered\ when
prompted\ \"Convert\ to\ file:\".

The\ other\ option\ \"Krigifier\ format\ ->\ gui\ format\"\ works\ in\ reverse,\ converting\ the
files\ intended\ to\ be\ read\ by\ the\ krigifier\ into\ gui-readable\ formats\ that\ may\ be
loaded\ using\ \"Load\".\ \ You\ will\ need\ to\ manually\ set\ the\ seed\ and\ stream\ in\ order
to\ reproduce\ functions.\ \ Further,\ due\ to\ the\ more\ varied\ nature\ of\ krigifier\ input,
this\ option\ is\ not\ as\ robust.\ \ Be\ warned."

   
    help_win $helptext
}

proc paramhelp {} {
    global WTitle
    global helptext
    set WTitle {Krigifier Parameters}

set helptext "p\tThe\ dimension\ of\ the\ space\ /\ The\ number\ of\ variables

n\tThe\ number\ of\ interpolated\ points.\ \ Larger\ values\ of\ n\ tend\ to
\tproduce\ functions\ with\ more\ minimizers.\ \ However,\ if\ n\ gets\ too\ large,
\tthe\ resulting\ functions\ can\ be\ fairly\ smooth.

tau^2\tMeasures\ the\ steepness\ of\ the\ underlying\ quadratic\ trend.\ \ The
\tratio\ of\ tau^2\ to\ sigma^2\ determines\ the\ strength\ of\ the\ noise\ over
\tthe\ trend.

beta0\tThe\ constant\ term\ of\ the\ underlying\ trend,\ beta0\ helps\ determine
\tthe\ likely\ height\ of\ the\ global\ minimum.

alpha\tDetermines\ the\ smoothness\ of\ the\ function.\ \ alpha\ =\ 1\ produces
\tfunctions\ which\ are\ everywhere\ differentiable,\ whereas\ alpha\ =\ 2\ produces
\tfunctions\ with\ sharp\ peaks\ at\ the\ interpolated\ points.

theta\tThe\ strength\ of\ the\ correlation.\ \ Large\ theta's\ produce\ functions\ with
\tmany\ sudden\ peaks\ and\ valleys.\ \ Small\ theta's\ result\ in\ more\ gradually\ changing
\tfunctions.

sigma^2\tThe\ variance\ of\ the\ noise.\ \ Larger\ values\ tend\ to\ produce\ greater
\tdeviations\ from\ the\ trend.

lower\t 
upper\tThese\ values\ set\ the\ boundaries\ of\ the\ design\ space.\ \ For\ example
\tlower\ =\ 0\ and\ upper\ =\ 1\ would\ be\ the\ unit\ cube." 
   
    help_win $helptext
}

proc versionhelp {} {
    global version
    global WTitle
    global helptext
    set WTitle {Krigifier Parameters}

set helptext "$version

The\ Krigifier\ Graphical\ Interface\ was\ written\ using\ tk/tcl,\ borrowing
from\ Bill\ Bynum's\ submit\ script.

The\ Krigifier\ itself\ is\ written\ in\ C++\ using\ an\ object-oriented\ approach.

All\ software\ is\ available\ from\ http://www.cs.wm.edu/~va/software

\tPermission\ to\ use,\ copy,\ modify,\ and\ distribute\ this\ software
\tfor\ any\ purpose\ without\ fee\ is\ hereby\ granted,\ provided\ that
\tthis\ entire\ notice\ is\ included\ in\ all\ copies\ of\ any\ software
\twhich\ is\ or\ includes\ a\ copy\ or\ modification\ of\ this\ software
\tand\ in\ all\ copies\ of\ the\ supporting\ documentation\ for\ such
\tsoftware.\ \ THIS\ SOFTWARE\ IS\ BEING\ PROVIDED\ \"AS\ IS\",\ WITHOUT
\tANY\ EXPRESS\ OR\ IMPLIED\ WARRANTY.\ \ IN\ PARTICULAR,\ THE\ AUTHOR
\tOFFERS\ NO\ REPRESENTATION\ OR\ WARRANTY\ OF\ ANY\ KIND
\tCONCERNING\ THE\ MERCHANTABILITY\ OF\ THIS\ SOFTWARE\ OR\ ITS
\tFITNESS\ FOR\ ANY\ PARTICULAR\ PURPOSE.


author:\ \ Anthony\ Padula\ \ adpadu@wm.edu"

 help_win $helptext
}


# The main stuff

wm title . "The Krigifier: A Pseudorandom Function Generator"

# Frame for the menu
frame .mf
pack .mf -side top -expand 1 -fill x
menubutton .mf.file -text "File" -menu .mf.file.m -width 20\
	-relief raised -anchor w
menu .mf.file.m -tearoff 0
pack .mf.file -side left -fill both -expand 1 -anchor w
.mf.file.m add checkbutton -label {Sound} -variable soundflag -onvalue 0 -offvalue 1
.mf.file.m add separator 
.mf.file.m add command -label {Save} -command {saveparams}
.mf.file.m add command -label {Load} -command {loadparams}
.mf.file.m add separator
.mf.file.m add command -label {Gui format -> krigifier format} -command {convert 0}
.mf.file.m add command -label {Krigifier format -> Gui format} -command {convert 1}
.mf.file.m add separator
.mf.file.m add command -label {Exit} -command {puts stdout "-1"; flush stdout; exit}

menubutton .mf.help -text "Help" -menu .mf.help.m -width 20\
	-relief raised -anchor w
menu .mf.help.m -tearoff 0
pack .mf.help -side left -fill both -expand 1 -anchor w
.mf.help.m add command -label {Help on Data Conversion} -command {converthelp}
.mf.help.m add command -label {Help on Parameters} -command {paramhelp}
.mf.help.m add command -label {About the Krigifier} -command {versionhelp}
# .mf.help.m.par -text "Help on Parameters" -menu .mf.help.m.par.m\
#	-width 20 -relief raised -anchor w -direction right
#menu .mf.help.m.par -tearoff 0
#pack .mf.help.m.par -side left -fill both -expand 1 -anchor w
#.mf.help.m.par add command -label {p} -command {puts stdout $p}

pack .mf.file -side left -in .mf
pack .mf.help -side left -in .mf

# Frame for holding buttons
frame .bf
pack  .bf -side bottom -expand 1 -fill x

# Frame to hold seed and stream
frame .sdf
pack .sdf -side bottom -expand 1 -fill x
make_entry_frame . sdf.seedf 0 {Seed (use -1 for random seed): } 30 {-1} 12
make_entry_frame . sdf.streamf 1 {Stream (0 - 255): } 20 {42} 12
pack .sdf.streamf -side right -expand 1 -fill x
pack .sdf.seedf -side left -expand 1 -fill x

# Frame to hold status
canvas .cs -width 600 -height 60 -bg white
pack .cs -side bottom

# Frame to hold sliders
frame .sf
pack .sf -side bottom -expand 1 -fill x

# Exit button
#button .bf.exit -text "Exit" -command {puts stdout "-1"; flush stdout; exit} -state active

# sound button
#button .bf.sound -text "Silent" -command {changesound}

# p button
radiobutton .bf.pbut1 -text "p = 1" -variable p -value 1 -indicatoron 1
radiobutton .bf.pbut2 -text "p = 2" -variable p -value 2 -indicatoron 1


# alpha button
radiobutton .bf.alphabut1 -text "alpha = 1" -variable alpha -value 1 -indicatoron 1
radiobutton .bf.alphabut2 -text "alpha = 2" -variable alpha -value 2 -indicatoron 1

# replot button
button .bf.replot -text "Plot" -command {replot}

# Save button
#button .bf.save -text "Save" -command {saveparams}

# Pack buttons into frame
#pack .bf.exit .bf.sound .bf.pbut .bf.alphabut .bf.save .bf.replot -side left -expand 1 -fill x
pack .bf.pbut1 .bf.pbut2 .bf.alphabut1 .bf.alphabut2 .bf.replot -side left -expand 1 -fill x

#put in text
.cs create text 10 15 -anchor w -fill red\
 -font deffont\
 -text $status -tag stattag

.cs create text 10 45 -anchor w -fill black\
 -font deffont\
 -text $fname -tag ftag

#.cs create text 10 75 -anchor w -fill black\
# -font "-adobe-helvetice-medium-o-normal--10"\
# -text $seedstring -tag stag

# n slider
scale .sf.n -from 0 -to 1000 -label "n" -length 600 -showvalue 1\
 -tickinterval 50 -variable n -width [expr $barwid] -orient horizontal

# tau2 slider
scale .sf.tau2 -from 0 -to 100 -label "tau^2" -length 600 -showvalue 1\
 -tickinterval 5 -variable tau2 -width [expr $barwid] -orient horizontal

# beta0 slider
scale .sf.beta0 -from -1000 -to 1000 -label "beta0" -length 600 -showvalue 1\
 -tickinterval 100 -variable beta0 -width [expr $barwid] -orient horizontal

# theta slider
scale .sf.theta -from 0 -to 1000 -label "theta" -length 600 -showvalue 1\
 -tickinterval 50 -variable theta -width [expr $barwid] -orient horizontal

# sigma2 slider
scale .sf.sigma2 -from 0 -to 1000 -label "sigma^2" -length 600 -showvalue 1\
 -tickinterval 50 -variable sigma2 -width [expr $barwid] -orient horizontal

# lower slider
scale .sf.lower -from -100 -to 100 -label "lower" -length 600 -showvalue 1\
 -tickinterval 10 -variable lower -width [expr $barwid] -orient horizontal\
 -command { changelower }

# upper slider
scale .sf.upper -from -100 -to 100 -label "upper" -length 600 -showvalue 1\
 -tickinterval 10 -variable upper -width [expr $barwid] -orient horizontal\
 -command { changeupper }

pack .sf.n .sf.tau2 .sf.beta0 .sf.theta .sf.sigma2 .sf.lower .sf.upper\
 -side top -expand 1 -fill x


