
breed [predators predator]
breed [preys prey]
breed [s1s s1]
breed [s0s s0]
breed [posteriors posterior]
breed [neurones neurone]
breed [spikes spike]
breed [ghosts ghost]
breed [bars bar]
breed [avgs avg]
breed [std-ups std-up]
breed [std-downs std-down]

s1s-own [ps1]
s0s-own [ps0]
posteriors-own [actual]
ghosts-own [actual]
std-ups-own [actual]
std-downs-own [actual]
avgs-own [actual]

globals [dt scale scale-avg marker marker-avg draw S resolution ghost-trail expect standard convergence-time]

to setup

  ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  __clear-all-and-reset-ticks
  resize-world -16 16 -18 16
  set convergence-time []
  reset

end

to reset

  clear-turtles

  set-patch-size 15
  ask patches [set pcolor white]

  set resolution 400

  create-s1s resolution [set shape "dot" set size .5 set heading 90 set color black
    setxy (-14 + who * 28 / count s1s) -18]

  create-s0s resolution [set shape "dot" set size .5 set color black
    setxy (-14 + (who - resolution) * 28 / resolution) -18]

    create-posteriors resolution [
      set shape "dot" set size .3 set color 0 setxy (-14 + (who - 2 * resolution) * 28 / resolution) 0
    ]


  ask turtles with [(breed = s1s or breed = s0s or breed = posteriors) and (xcor < -12.3)] [die]

  if prior? = "normal" [
    let mu .25
    let sigma .1
    let A (1 / (sigma * sqrt(2 * pi)))
    ask posteriors [
      set actual (A * e ^ (-1 * ((([xcor] of self + 14) / 28) - (mu + 14) / 28) ^ 2 / (2 * sigma ^ 2)) )
    ]
  ]

  if prior? = "dirac-delta" [
    ask min-one-of posteriors [abs xcor] [set actual resolution]
  ]

  if prior? = "uniform" [
    ask posteriors [set actual resolution / count s1s]
  ]



  set dt .001

  ask s0s [let I (g * (([xcor] of self + 14) / 56) ^ (-2) + n) * dt
    set ps0 e ^ (-1 * I)
    set ycor (ycor + 4 * e ^ (-1 * I))]

  ask s1s [let I (g * (([xcor] of self + 14) / 56) ^ (-2) + n) * dt
    set ps1 1 - e ^ (-1 * I)
    set ycor (ycor + 4 * (1 - e ^ (-1 * I)))
    ]

  create-predators 1 [setxy -14 15 set size 2.5 set shape "circle" set color red]
  create-neurones 1 [setxy -14 15 set size 1 set shape "circle" set color black]
  create-preys 1 [setxy (-14 + D * 28 * 2) 15 set size 1 set shape "circle" set color blue]

  ask patch -13 -14 [set plabel-color black set plabel "Pr(S=1|D)"]
  ask patch -13 -18 [set plabel-color black set plabel "Pr(S=0|D)"]

  set scale 14 / (2 * floor ([actual] of one-of posteriors with-max [actual]))
  set marker round (14 / scale)
  ask patch -15 12 [set plabel-color black set plabel marker]
  ask patch -15 -13 [set plabel-color black set plabel 0]
  ask patch 13 -13 [set plabel-color black set plabel "red = posterior mean, blue = +/- 3 posterior std. devs."]
  ask posteriors [set ycor (0 + scale * actual)]
  if prior? = "dirac-delta" [ask posteriors with [xcor = 0] [set shape "line" set heading 0 set ycor ycor - 3.5 set size scale * actual]]

  ;get probabilities of spiking
  let signal (1 / d ^ 2)
  let intensity (g * signal + n)
  set draw (1 - e ^ (- intensity * dt))

  create-bars 1 [set color black set shape "line" set size 32 set heading 90 setxy 0 13.5]
  create-bars 1 [set color black set shape "line" set size 32 set heading 90 setxy 0 -.5]
  create-bars 1 [set color black set shape "line" set size 32 set heading 90 setxy 0 -1.5]
  create-bars 1 [set color black set shape "line" set size 32 set heading 90 setxy 0 -13.5]

end

to go

  ask spikes [set xcor xcor - .5]

    set S random-float 1

    ifelse S < draw [
    ask neurones [set color white]
    display
    wait .05
    ask neurones [set color black]
    create-spikes 1 [set shape "line" set size .7 setxy 14 -1.1 set color black set heading 0]
    display
    ]
    [ask neurones [set color black]
    ]

  ask spikes with [xcor < -15] [die]

  if (prior? = "normal" or prior? = "uniform") [update-posterior]

  tick

end

to update-posterior

  ask ghosts [set color color + .9]
  ask ghosts with [color > 9.9] [die]

  ask posteriors [hatch-ghosts 1 [
      set size .3 set color 1 set shape "dot"]]

  ;update the posterior
  ifelse S < draw
    [ask posteriors [set actual (actual * item 0 [ps1] of s1s with [xcor = [xcor] of myself])]]
    [ask posteriors [set actual (actual * item 0 [ps0] of s0s with [xcor = [xcor] of myself])]]

  ;normalise new posterior
  let normalise (sum [actual] of posteriors)
  ask posteriors [set actual (actual * resolution / normalise)]

  plot-smart

end

to plot-smart

  ;does the plot need to be rescaled?
  let top max [ycor] of turtles with [breed = posteriors or breed = ghosts]

  ifelse (top > 11 or top < 6)
    [set scale 13 / (2 * floor ([actual] of one-of turtles with [breed = posteriors or breed = ghosts] with-max [actual]))
     set marker round (13 / scale)
     if (((marker < 50) and ((prior? = "uniform") or (prior? = "normal"))) or (marker < 500 and prior? = "dirac-delta")) [
     ask patch -15 12 [set plabel-color black set plabel marker]
     ask ghosts [set ycor (0 + scale * actual)]]
     ask posteriors [set ycor (0 + scale * actual)]
    ]

    [ask ghosts [set ycor (0 + scale * actual)]
     ask posteriors [set ycor (0 + scale * actual)]]

  ;get mean and spread of posterior and plot it
  ask avgs [set xcor (xcor - .5)]
  ask std-ups [set xcor (xcor - .5)]
  ask std-downs [set xcor (xcor - .5)]

  set expect sum [actual * (xcor + 14) / 56] of posteriors / resolution
  set standard (((sum [actual * ((xcor + 14) / 56) ^ 2] of posteriors / resolution) - expect ^ 2) ^ (1 / 2)) * 3

  create-avgs 1 [set color red set shape "dot" set size .5 set actual expect setxy 14 actual]
  create-std-ups 1 [set color blue set shape "dot" set size .3 set actual (expect + standard) setxy 14 actual]
  ifelse (expect - standard < 0)
    [create-std-downs 1 [set color blue set shape "dot" set size .3 set actual 0 setxy 14 actual]]
    [create-std-downs 1 [set color blue set shape "dot" set size .3 set actual (expect - standard) setxy 14 actual]]

  let upper max [ycor] of std-ups
  let lower min [ycor] of std-downs

  ifelse (upper > -4 or upper < -8)
    [set scale-avg 11 / ([actual] of one-of std-ups with-max [actual])
     set marker-avg (round (1000 * 11 / scale-avg)) / 1000
     ask patch -15 -2 [set plabel-color black set plabel marker-avg]
     ask avgs [set ycor (-13 + scale-avg * actual)]
     ask std-ups [set ycor (-13 + scale-avg * actual)]
     ask std-downs [set ycor (-13 + scale-avg * actual)]
    ]
    [ask avgs [set ycor (-13 + scale-avg * actual)]
     ask std-ups [set ycor (-13 + scale-avg * actual)]
     ask std-downs [set ycor (-13 + scale-avg * actual)]]

  ask avgs with [xcor < -15] [die]
  ask std-ups with [xcor < -15] [die]
  ask std-downs with [xcor < -15] [die]





end
@#$#@#$#@
GRAPHICS-WINDOW
366
10
871
566
16
-1
15.0
1
10
1
1
1
0
1
1
1
-16
16
-18
16
1
1
1
ticks
30.0

BUTTON
146
352
212
385
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
129
123
221
156
D
D
.01
.3
0.1
.01
1
NIL
HORIZONTAL

SLIDER
129
164
221
197
g
g
1
15
6
1
1
NIL
HORIZONTAL

SLIDER
129
204
221
237
n
n
0
100
100
1
1
NIL
HORIZONTAL

CHOOSER
108
249
246
294
prior?
prior?
"normal" "uniform" "dirac-delta"
1

BUTTON
145
481
213
514
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
64
52
313
112
1)  Set distance between predator and prey, parameters of the sensor, and initial prior
16
0.0
1

TEXTBOX
124
322
274
342
2)  Click setup
16
0.0
1

TEXTBOX
81
423
296
464
3)  Click go; click once to run, click again to pause
16
0.0
1

TEXTBOX
6
20
358
51
Note:  see Information tab for more details and observations
11
0.0
1

TEXTBOX
67
535
310
564
A speed slider can be found at the top of the screen to adjust the speed of the simulation
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This model illustrates how individual spikes and non-spikes can be interpreted as   measurements about some state of the world, and how we can perform Bayesian inference on those measurements to infer about the world.  In this model, the state of the world considered is the distance D from a predator (red) to a prey (blue).  Given a probabilistic model of spiking given D, Pr(S|D), Bayes' rule states how to calculate Pr(D|S) given the outputs of the neurone.

## HOW IT WORKS

The predator has a neurone that fires a spike on a time step with probability Pr(S|D), as outlined in the manuscript.  Given this conditional probability, which is plotted at the bottom of the screen, and a prior distribution Pr(D), we can infer the prey location Pr(D|S) by Bayes' rule.  By letting the posterior distribution Pr(D|S) become the prior on the next time step, we can perform Bayes' rule recursively in time.

## HOW TO USE IT

To start, set the parameters of the model to your desired values.  Choose whether the prior is a normal, uniform, or Dirac-delta distribution.  Click the `setup' button, and you're ready to run simulations.

Next, click the `go' button.  Click the button once to start running a simulation, and click again to pause it.  The model is run at a sufficiently slow speed so that the viewer can see how individual spikes (and non-spikes) can be used to perform Bayesian inference about D.  Every time step is one millisecond in duration.  When the predator's neurone spikes in a time step, the predator's neurone flashes white and makes a pulse of sound.  By Bayes' rule, the posterior Pr(D|S) shifts to the left.  When the neurone doesn't spike in a time step, the posterior decays to the right.  The observed spike train is streamed below.  The sensor's conditional probabilities of spiking are in the lowest plot.

## THINGS TO NOTICE

Pr(D|S) converges quickly to the observed prey location when D is small, often in less than 500 time steps (one-half of one second in real-time), and even when the gain/noise ratio for the sensor is extremely small.  When D is small (say 0.05 or so), then Pr(D|S) quickly converges on D = d, in as little as 20 or 30 timesteps, despite low gain-to-noise ratios.

The convergence of the posterior distribution to the observed prey location is insensitive to the prior distribution of the prey, so long as the prior has reasonably heavy tails.  To verify this, run the model for a uniform, normal, and dirac-delta prior by changing the chooser labeled 'prior?'  The posterior Pr(D|S) will converge on the observed prey location D = d, so long as the prior is not inconsistent with the observed prey location.

## THINGS TO TRY

Investigate how the gain-to-noise ratio affects the speed of convergence of Pr(D|S) to the observed value of D = d.

Investigate how the observed value of D = d affects the speed of convergence of Pr(D|S).  The larger D = d is, the longer it takes Pr(D|S) to converge to D = d.

## EXTENDING THE MODEL

Spatial distributions of these sensors can provide information on the location of other organisms in higher dimensions.  This is consistent with the notion of "skin brains" (Holland 2003).

To deal with the more general case of prey moving, we need to introduce the conditional probability Pr(X(t)|X(t-1)) (a dynamic prior).  Then, the neurone becomes a rudimentary particle filter.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
