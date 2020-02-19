# Octave example

```bash
cd Stalefish/ucr
octave # Run the GNU Octave program
```
These are example interactive commands in octave. They can be scripted
if necessary.

```
octave:1> load V_Id2_1.h5
octave:2> who
Variables in the current scope:

Frame000  Frame007  Frame014  Frame021  Frame028  Frame035  Frame042
Frame001  Frame008  Frame015  Frame022  Frame029  Frame036  Frame043
Frame002  Frame009  Frame016  Frame023  Frame030  Frame037  Frame044
Frame003  Frame010  Frame017  Frame024  Frame031  Frame038  Frame045
Frame004  Frame011  Frame018  Frame025  Frame032  Frame039  Frame046
Frame005  Frame012  Frame019  Frame026  Frame033  Frame040  ans
Frame006  Frame013  Frame020  Frame027  Frame034  Frame041  nframes

octave:3> Frame000
Frame000 =

  1x1 struct array containing the fields:

    box0
    box1
... (snip) ...
    box98
    box99
    class
    fitted
    fitted_offset
    fitted_rotated
    means
    means_autoscaled
    sbox_centers
    sbox_linear_distance
    sboxes

octave:4> Frame000.class
ans =

  1x1 struct array containing the fields:

    P
    PP000
    PP001
    PP002
    PP_n
    binA
    binB
    filename
    flags
    idx
    layer_x
    nBinsTarg
    pixels_per_mm
    polyOrder
    pp_idx
    thickness

octave:5> Frame000.class.nBinsTarg
ans = 100
octave:6> Frame000.class.pixels_per_mm
ans =    233
octave:7> Frame000.class.PP000
ans =

  507  462  454  473  506
  617  554  487  427  378

octave:8>
```
