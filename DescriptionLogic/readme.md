# Breast Cancer Treatment Staging Using Description Logic

## Overview

The [`breast-cancer.owl`](breast-cancer.owl) file contains the ontology that models the domain of
breast cancer and uses restriction classes to classify example diagnoses into
treatment categories. The file was generated using Protégé 5.1.0 and tested using
the HermiT reasoner.

## Domain Model

The ontology contains roughly 200 classes (see [full list](https://github.com/pmanko/BIME550_Project4/blob/master/DescriptionLogic/class-hierarchy.txt)), which
represent numerous features related to the patient, the diagnosis and treatment
options.

## Inference

The ontology can be used to answer two main competency questions, described below.

### Surgery

The first question is regarding which type of surgery a patient should get after
diagnosis. We have implemented both ways that this decision can be made. First,
if the cancer has been staged, we can use this as follows:

```
hasStage some (StageIIA or StageIIB or StageIIIB or StageIIIC or StageIV)
```

The second method implemented uses the outcome of preoperative systemic therapy to
classify the test examples (see NCCN decision tree number BINV-12).

We have also implemented the staging restriction classes as in the previous
assignment, but specific to breast cancer.

### Adjuvant Therapy

Adjuvant therapy is one or more treatments given after surgery. The NCCN has
many complicated decision trees to make this decision. Our ontology implements a
simplified version of these rules (NCCN decision tree numbers BINV-4 to BINV-9).

An example of our simplified restriction class logic is:

```
(
  (hasHistology some (DuctalHistology or LobularHistology or MetaplasticHistology or MixedHistology or MucinousHistology or TubularHistology))
    and
  ((hasERStatus some HRPositive) or (hasPRStatus some HRPositive))
    and
  (((hasHER2Status some HER2Positive) and ((hasDistalMetastasis some M1) or ((hasDistalMetastasis some M0) and (hasPrimaryClassification some (TumorSizeGreaterThan5mmLessThan10mm or TumorSizeGreaterThan10mm)))))
     or
   ((hasHER2Status some HER2Negative) and (hasDistalMetastasis some M0) and ((hasPrimaryClassification some (TumorSizeGreaterThan5mmLessThan10mm or TumorSizeGreaterThan10mm)) and (hasRecurrance some (NoGeneticRecurrance or IntermediateGeneticRecurrance))))
  )
)
```

## Testing

To test the restriction classes, load the file into Protégé and run the
reasoner. Then, navigate to the `Diagnosis > TestDiagnoses` parent class and
browse through its children. In each case there should be one or more
inferred subclasses that match the name of the test diagnosis.
