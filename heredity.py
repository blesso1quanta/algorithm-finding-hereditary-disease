import csv
import itertools
import sys

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    jointprobability = {}
    for i in people:
        if i in two_genes:
            if people[i]["mother"] == None and people[i]["father"] == None:
                if i not in have_trait:
                    probability1 = PROBS["gene"][2]
                    probability2 = PROBS["trait"][2][False]
                    jointprobability[i] = probability2 * probability1
                else:
                    probability1 = PROBS["gene"][2]
                    probability2 = PROBS["trait"][2][True]
                    jointprobability[i] = probability2 * probability1
            else:
                # the gene from mother and from his father
                if people[i]["father"] in one_gene:
                    probability3 = 0.505
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.99
                    else:
                        probability4 = 0.01
                    prejointprobality = probability3 * probability4
                elif people[i]["father"] in two_genes:
                    probability3 = 0.99
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.99
                    else:
                        probability4 = 0.01
                    prejointprobality = probability3 * probability4

                else:
                    probability3 = 0.01
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.99
                    else:
                        probability4 = 0.01
                    prejointprobality = probability3 * probability4
                if i not in have_trait:
                    probability7 = PROBS["trait"][2][False]
                    jointprobability[i] = prejointprobality * probability7
                else:
                    probability8 = PROBS["trait"][2][True]
                    jointprobability[i] = prejointprobality * probability8


        elif i in one_gene:
            if people[i]["mother"]==None and people[i]["father"]==None:
                if i not in have_trait:
                    probability1 = PROBS["gene"][1]
                    probability2 = PROBS["trait"][1][False]
                    jointprobability[i] = probability2 * probability1
                else:
                    probability1 = PROBS["gene"][1]
                    probability2 = PROBS["trait"][1][True]
                    jointprobability[i] = probability2 * probability1
            else:
                #the gene from mother and not from his father
                if people[i]["father"] in one_gene:
                    probability3 = 0.505
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.99
                    else:
                        probability4 = 0.01
                    prejointprobality1 = probability3 * probability4
                elif people[i]["father"] in two_genes:
                    probability3 = 0.01
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.99
                    else:
                        probability4 = 0.01
                    prejointprobality1 = probability3 * probability4

                else:
                    probability3 = 0.99
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.99
                    else:
                        probability4 = 0.01
                    prejointprobality1 = probability3 * probability4

                #he gets the gene from his father and not from his mother
                if people[i]["mother"] in one_gene:
                    probability5 = 0.505
                    if people[i]["father"] in one_gene:
                        probability6 = 0.505
                    elif people[i]["father"] in two_genes:
                        probability6 = 0.99
                    else:
                        probability6 = 0.01
                    prejointprobality2 = probability5 * probability6
                elif people[i]["mother"] in two_genes:
                    probability5 = 0.01
                    if people[i]["father"] in one_gene:
                        probability6 = 0.505
                    elif people[i]["father"] in two_genes:
                        probability6 = 0.99
                    else:
                        probability6 = 0.01
                    prejointprobality2 = probability5 * probability6
                else:
                    probability5 = 0.99
                    if people[i]["father"] in one_gene:
                        probability6 = 0.505
                    elif people[i]["father"] in two_genes:
                        probability6 = 0.99
                    else:
                        probability6 = 0.01
                    prejointprobality2 = probability5 * probability6
                prejointprobality = prejointprobality2 + prejointprobality1


                if i not in have_trait:
                    probability7 = PROBS["trait"][1][False]
                    jointprobability[i] = prejointprobality * probability7
                else:
                    probability8 = PROBS["trait"][1][True]
                    jointprobability[i] = prejointprobality * probability8
        else:
            if people[i]["mother"] == None and people[i]["father"] == None:
                if i not in have_trait:
                    probability1 = PROBS["gene"][0]
                    probability2 = PROBS["trait"][0][False]
                    jointprobability[i] = probability2 * probability1
                else:
                    probability1 = PROBS["gene"][0]
                    probability2 = PROBS["trait"][0][True]
                    jointprobability[i] = probability2 * probability1
            else:
                # the gene not from mother and not from his father
                if people[i]["father"] in one_gene:
                    probability3 = 0.505
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.01
                    else:
                        probability4 = 0.99
                    prejointprobality = probability3 * probability4
                elif people[i]["father"] in two_genes:
                    probability3 = 0.01
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.01
                    else:
                        probability4 = 0.99
                    prejointprobality = probability3 * probability4

                else:
                    probability3 = 0.99
                    if people[i]["mother"] in one_gene:
                        probability4 = 0.505
                    elif people[i]["mother"] in two_genes:
                        probability4 = 0.01
                    else:
                        probability4 = 0.99
                    prejointprobality = probability3 * probability4
                if i not in have_trait:
                    probability7 = PROBS["trait"][0][False]
                    jointprobability[i] = prejointprobality * probability7
                else:
                    probability8 = PROBS["trait"][0][True]
                    jointprobability[i] = prejointprobality * probability8
    totalprobability = 1
    for i in jointprobability:
        totalprobability = totalprobability * jointprobability[i]


    return totalprobability

def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for i in probabilities:
        if i in two_genes:
            probabilities[i]["gene"][2] = p
            if i in have_trait:
                probabilities[i]["trait"][True] = p
            else:
                probabilities[i]["trait"][False] = p
        elif i in one_gene:
            probabilities[i]["gene"][1] = p
            if i in have_trait:
                probabilities[i]["trait"][True] = p
            else:
                probabilities[i]["trait"][False] = p
        else:
            probabilities[i]["gene"][0] = p
            if i in have_trait:
                probabilities[i]["trait"][True] = p
            else:
                probabilities[i]["trait"][False] = p





def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for i in probabilities:
        t = probabilities[i]["trait"][True]
        f = probabilities[i]["trait"][False]
        if t>=f:
            if f!=0:
                nfactor = round(t / f)
                ffactor = 1 / (nfactor + 1)
                tfactor = ffactor * nfactor
                probabilities[i]["trait"][True] = tfactor
                probabilities[i]["trait"][False] = ffactor
            else:
                probabilities[i]["trait"][True] = 1.0


        else:
            if t!=0:
                nfactor = round(f / t)
                tfactor = 1 / (nfactor + 1)
                ffactor = tfactor * nfactor
                probabilities[i]["trait"][True] = tfactor
                probabilities[i]["trait"][False] = ffactor
            else:
                probabilities[i]["trait"][False] = 1.0






if __name__ == "__main__":
    main()
