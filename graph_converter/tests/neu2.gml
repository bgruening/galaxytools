graph [
  directed 1
  node [
    id 1
    label "A"
    weiteresAttribut "42"
  ]
  node [
    id 3
    label "C"
    weiteresAttribut "44"
  ]
  node [
    id 2
    label "B"
    weiteresAttribut "43"
  ]
  edge [
    source 1
    target 2
    weiteresAttribut "44"
    id "3"
    label "Kante AB"
  ]
  edge [
    source 3
    target 1
    weiteresAttribut "44"
    id "3"
    label "Kante CA"
  ]
  edge [
    source 2
    target 3
    weiteresAttribut "44"
    id "3"
    label "Kante BC"
  ]
]
