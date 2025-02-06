(define (problem minimal-bug-case)
   (:domain miconic)
   (:objects p0 p1 p2
             f0 f1 f2 f3 f4)

(:init
(passenger p0)
(passenger p1)
(passenger p2)

(floor f0)
(floor f1)
(floor f2)
(floor f3)
(floor f4)

(above f0 f1)
(above f0 f2)
(above f0 f3)
(above f0 f4)

(above f1 f2)
(above f1 f3)
(above f1 f4)

(above f2 f3)
(above f2 f4)

(above f3 f4)

(origin p0 f1)
(destin p0 f3)

(origin p1 f4)
(destin p1 f2)

(origin p2 f0)
(destin p2 f4)

(lift-at f0)
)

(:goal (and
(served p0)
(served p1)
(served p2)
))
)
