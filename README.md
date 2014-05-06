ARGBEAST
========

ARGBEAST is a [BEAST 2](http://beast2.org) package which facilitates
inference of a (restricted kind of) ancestral recombination graph
(ARG) and related parameters from a sequence alignment.  It is
inspired by
[Didelot et al.'s 2010 Genetics paper](http://www.genetics.org/content/186/4/1435)
and is currently under heavy development.  Check back soon!

[![Build Status](https://travis-ci.org/CompEvol/ARGBEAST.svg?branch=master)](https://travis-ci.org/CompEvol/ARGBEAST)


Road map for ARGBEAST development
---------------------------------

- [x] ARG state node
- [x] ARG likelihood
- [x] ARG density under model
- [x] ARG simulator

- [x] Operators

  - [x] Add/remove
  - [x] ARG scaler
  - [x] Clonal frame operators
  - [x] Converted region slide
  - [x] Converted region swap
  - [x] Recombinant edge slide
  - [x] Recombinant edge hop
  - [ ] Recombination merge/split *IN PROGRESS*
  - [ ] Clonal frame / recomb topology switch

- [x] Sequence alignment simulator
- [ ] Inference from simulated data *IN PROGRESS*
- [ ] Inference from real data
