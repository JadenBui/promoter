package qut.parallel;

import edu.au.jacobi.pattern.Match;
import qut.GenbankRecord;
import qut.Gene;
import qut.NucleotideSequence;


/**
 * A data container for the streaming work.
 *
 * @author Jordi Smit on 11-10-2018.
 */

public class GeneComparisionTask {
    private final NucleotideSequence nucleotideSequence;
    private final Gene gene;
    private final Gene referenceGene;
    private Match prediction;

    public GeneComparisionTask(NucleotideSequence nucleotideSequence, Gene referenceGene, Gene gene) {
        this.nucleotideSequence = nucleotideSequence;
        this.referenceGene = referenceGene;
        this.gene = gene;
        this.prediction = null;
    }

    public Match getPrediction() {
        return prediction;
    }

    public void setPrediction(Match prediction) {
        this.prediction = prediction;
    }

    public NucleotideSequence getNucleotideSequence() {
        return nucleotideSequence;
    }

    public Gene getGene() {
        return gene;
    }

    public Gene getReferenceGene() {
        return referenceGene;
    }
}
