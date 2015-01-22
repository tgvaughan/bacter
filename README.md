BACTER
======

BACTER is a [BEAST 2](http://beast2.org) package which facilitates
inference of a (restricted kind of) ancestral recombination graph
(ARG) and related parameters from a sequence alignment.  It is
inspired by
[Didelot et al.'s 2010 Genetics paper](http://www.genetics.org/content/186/4/1435)
and is currently under heavy development.  Check back soon!

[![Build Status](https://travis-ci.org/CompEvol/BACTER.svg?branch=master)](https://travis-ci.org/CompEvol/BACTER)

Development news
----------------

An initial version of BACTER is functional and reliably infers ARGs from
data simulated under a restricted ClonalOrigin-style model which assumes
that each site is affected by at most one conversion event.  However, we
are finding that this model is too restrictive for application to real data.

Therefore, we are currently in the process of lifting this restriction.


New development road map
------------------------

Phase 1:
- [x] Update ARG state node
- [x] Update existing ARG density, likelihood and simulation code to
       work with the new state
- [x] Update operators to work with new state
- [x] Reorganize packages so as to keep code specific to restricted model separate from new code.

Phase 2:
- [ ] ARG simulator for new model (allow for model misspec analysis) *IN PROGRESS*
- [ ] ARG density under new model
- [ ] ARG operators for new model


Old road map (legacy)
---------------------

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
  - [x] Recombination merge/split
  - [x] Clonal frame / recomb topology switch

- [x] Sequence alignment simulator
- [ ] Inference from data simulated under approximate model *IN PROGRESS*
- [ ] Inference from data simulated under full model
- [ ] Inference from real data *IN PROGRESS*
