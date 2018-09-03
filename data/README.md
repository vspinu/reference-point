
Dataset and code available at https://github.com/vspinu/reference-point

dataset.csv
===========

The dataset is organized in a 'long' format, with each row describing one
question for one subject.

  - IdSubject       : arbitrary number identifying subjects
  - IdQuestion      : question number from 1 to 70
  - x1a x2a x3a x4a : outcomes of prospect A
  - p1a p2a p3a p4a : probabilities of prospect A
  - x1b x2b x3b x4b : outcomes of prospect B
  - p1b p2b p3b p4b : probabilities of prospect B
  - Rank            : rank of the question in the experiment for the subject (randomized between subjects)
  - UpDown          : 1 if prospect A was displayed at the top of the screen, 2 if at the bottom
  - Preference      : 1 if prospect A was preferred, 2 if B was preferred
  - TimeChoice      : date and time of the choice


questions.csv
=============

  - IdQuestion      : question number from 1 to 70
  - x1a x2a x3a x4a : outcomes of prospect A
  - p1a p2a p3a p4a : probabilities of prospect A
  - x1b x2b x3b x4b : outcomes of prospect B
  - p1b p2b p3b p4b : probabilities of prospect B

subjects.csv
============

  - IdSubject: arbitrary number identifying subjects
  - Age
  - Gender (sometimes not reported)
