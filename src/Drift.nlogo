extensions [sound]

breed [predators predator]
breed [preys prey]
breed [s1s s1]
breed [s0s s0]
breed [posteriors posterior]
breed [neurones neurone]
breed [spikes spike]
breed [reserves reserve]
breed [ghosts ghost]

s1s-own [ps1]
s0s-own [ps0]
posteriors-own [actual]
reserves-own [actual]
ghosts-own [actual]

globals [dt scale marker draw S resolution d drift g n]

to setup 
  
  ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  __clear-all-and-reset-ticks
  resize-world -16 16 -16 16
  set-patch-size 15
  ask patches with [pycor >= 12] [set pcolor white]
  
  ask patches with [(pycor <= 7) and (pycor >= -7)] [set pcolor white]
  ask patches with [pycor = 10 or pycor = 9] [set pcolor white]
  
  ask patches with [pycor <= -9] [set pcolor white]
  
  set resolution 200
  
  create-s1s resolution [set shape "dot" set size .5 set heading 90 set color black
    setxy (-14 + who * 28 / resolution) -15]
  
  create-s0s resolution [set shape "dot" set size .5 set color black
    setxy (-14 + (who - resolution) * 28 / resolution) -15]
  
  create-posteriors resolution [
      set shape "dot" set size .3 set color 0 setxy (-14 + (who - 2 * resolution) * 28 / resolution) -6
  ]
  
  ask turtles with [(breed = s1s or breed = s0s or breed = posteriors) and (xcor < -13.44)] [die]
  
    let mu .25
    let sigma .1
    let A (1 / (sigma * sqrt(2 * pi)))
    ask posteriors [set actual resolution / count s1s]
    
  create-reserves resolution [
    set shape "dot" set size .3 set color white
    set xcor ([xcor] of max-one-of s1s [xcor]) 
    set actual [actual] of one-of posteriors
  ]
    
  
  set drift 28 / resolution
  set dt .001
  set d 1
  set g 20
  set n 100
  
  ask s0s [let I (g * (([xcor] of self + 14) / 14) ^ (-2) + n) * dt
    set ps0 e ^ (-1 * I)
    set ycor (ycor + 4 * e ^ (-1 * I))]
  
  ask s1s [let I (g * (([xcor] of self + 14) / 14) ^ (-2) + n) * dt
    set ps1 1 - e ^ (-1 * I)
    set ycor (ycor + 4 * (1 - e ^ (-1 * I)))
    ]
  
  create-predators 1 [setxy -14 14 set size 3 set shape "circle" set color red]
  create-neurones 1 [setxy -14 14 set size 1 set shape "circle" set color black]
  create-preys 1 [setxy ([xcor] of max-one-of posteriors [xcor]) 14 set size 1 set shape "circle" set color blue]
  
  ;shave off pesky decimals
  ask turtles [set xcor (round (xcor * 100) / 100)]
  
  ask patch 11 -16 [set plabel-color black set plabel "Pr(S=1|D)"]
  ask patch 11 -10 [set plabel-color black set plabel "Pr(S=0|D)"]
  
  set scale 14 / (2 * floor ([actual] of one-of posteriors with-max [actual]))
  set marker round (14 / scale)
  ask patch -15 7 [set plabel-color black set plabel marker]
  ask turtles with [breed = posteriors or breed = reserves] [set ycor (-6 + scale * actual)]

  
end

to go
  
  if d < .05 [stop]
  move-prey
  see-a-spike?
  drift-posterior
  update-posterior
  plot-smart
  if [xcor] of one-of preys < -9 [wait .1]
  tick
  
end

to move-prey
  
  ask preys [setxy (xcor - drift) 14]
  
end

to see-a-spike?
  
  ;get probabilities of spiking
  set d ([xcor] of one-of preys - [xcor] of one-of predators) / 14
  let signal (1 / d ^ 2)
  let intensity (g * signal + n)
  set draw (1 - e ^ (- intensity * dt))
  
  ask spikes [set xcor xcor + .3]
  
  set S random-float 1
  
    ifelse S < draw [
    sound:play-drum "ACOUSTIC SNARE" 64
    ask neurones [set color white]
    create-spikes 1 [set shape "line" set size 1 setxy -14 9 set color black set heading 0]
    display
    ]
    [ask neurones [set color black]]
    
  ask spikes with [xcor > 14] [die]
  
end

to drift-posterior
 
  
  ask posteriors [set xcor (xcor - drift)]
  ask one-of reserves with-min [who] [
    set color black set breed posteriors set shape "dot"
    set xcor (xcor - drift)
  ]
  
  ask posteriors with [xcor < -13.44] [die]
  
  ;shave off pesky decimals
  ask turtles [set xcor (round (xcor * 100) / 100)]
  
end

to update-posterior
  
  ;update the newest posterior
  ifelse S < draw 
    [ask posteriors [set actual (actual * item 0 [ps1] of s1s with [xcor = [xcor] of myself])]
     ask reserves [set actual (actual * item 0 [ps1] of s1s with [xcor = [xcor] of myself])]
    ]
    [ask posteriors [set actual (actual * item 0 [ps0] of s0s with [xcor = [xcor] of myself])]
     ask reserves [set actual (actual * item 0 [ps0] of s0s with [xcor = [xcor] of myself])]
    ]
    
  ;normalise new posterior
  let normalise (sum [actual] of posteriors)
  ask posteriors [set actual (actual * resolution / normalise)]
  ask reserves [set actual (actual * resolution / normalise)]
  
  ask ghosts [set color color + .9]
  ask ghosts with [color > 9.9] [die]
  
  ask posteriors [hatch-ghosts 1 [
      set size .3 set color 1 set shape "dot"]]
  
end

to plot-smart
  
  ;does the plot need to be rescaled?
  let top max [ycor] of turtles with [breed = posteriors or breed = ghosts]
  
  ifelse ((top > 6) or (top < 0))
    [set scale 14 / (2.5 * floor ([actual] of one-of turtles with [breed = posteriors or breed = ghosts] with-max [actual]))
     set marker round (14 / scale)
     ask patch -15 7 [set plabel-color black set plabel marker]
     ask posteriors [set ycor (-6 + scale * actual)]
     ask reserves [set ycor (-6 + scale * actual)]
     ask ghosts [set ycor (-6 + scale * actual)]
    ]
    
    [ask posteriors [set ycor (-6 + scale * actual)]
     ask reserves [set ycor (-6 + scale * actual)]
     ask ghosts [set ycor (-6 + scale * actual)]
    ]
    
    ask patch -3 -7 [
      ifelse [xcor] of one-of preys > 0 
      [set plabel-color black set plabel "prey is almost certainly not very close"]
      [set plabel ""]
    ]
     
    ask patch 0 -7 [
      ifelse [xcor] of one-of preys < 0 and [xcor] of one-of preys > -8
      [set plabel-color black set plabel "as prey gets closer, watch posterior follow it in"]
      [set plabel ""]
    ]
    
    ask patch 0 -7 [
      if [xcor] of one-of preys < -13
      [set plabel-color black set plabel "now prey is almost certainly very close"]
    ]
    
    ask patch 14 -5 [
      if [xcor] of one-of preys < -13
      [set plabel-color black set plabel "<- fading trail indicates posterior tracking the prey when it was very close"]
    ]
      

  
end
@#$#@#$#@
GRAPHICS-WINDOW
346
26
851
552
16
16
15.0
1
10
1
1
1
0
0
1
1
-16
16
-16
16
1
1
1
ticks
30.0

BUTTON
135
114
201
147
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

BUTTON
135
213
203
246
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
107
86
257
106
1)  Click setup
16
0.0
1

TEXTBOX
70
160
285
201
2)  Click go; click once to run, click again to pause
16
0.0
1

TEXTBOX
73
35
303
53
Note:  see Information tab for more details
11
0.0
1

TEXTBOX
51
254
294
283
A speed slider can be found at the top of the screen to adjust the speed of the simulation
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This model illustrates how sequential Bayesian inference can be applied to models where the state of the world being measured is dynamic.

## HOW IT WORKS

The predator has a neurone that fires a spike on a time step with probability Pr(S|D), as outlined in the manuscript.  Given this conditional probability, which is plotted at the bottom of the screen, and a prior distribution Pr(D), we can infer the prey location Pr(D|S) by Bayes' rule.  By letting the posterior distribution Pr(D|S) become the prior on the next time step, we can perform Bayes' rule recursively in time.

However, the prey deterministically drifts closer to the predator after every time step.  When the state of the world is dynamic, we need to calculate the posterior using a two-step recursion as outlined in the text.  After a time step, the posterior shifts with the prey and is then updated given the spiking output of the neurone.  The process then repeats.

## HOW TO USE IT

Just click 'setup' and 'go.'

## THINGS TO NOTICE

The notes that appear during the model's running time point out key observations.
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
NetLogo 5.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 1.0 0.0
0.0 1 1.0 0.0
0.2 0 1.0 0.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
