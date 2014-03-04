graph [
	comment "Das ist ein Beispielgraph."
	directed 1
	id 42
	label "Graph"
	node [
		id 1
		label "A"
		weiteresAttribut 42
	]
	node [
		id 2
		label "B"
		weiteresAttribut 43
	]
	node [
		id 3
		label "C"
		weiteresAttribut 44
	]
	edge [
		source 1
		target 2
		label "Kante AB"
	]
	edge [
		source 2
		target 3
		label "Kante BC"
	]
	edge [
		source 3
		target 1
		label "Kante CA"
	]
]
