package qut.parallel;

import edu.au.jacobi.pattern.Match;

public class Prediction {
    private final Match predictionResult;
    private final String referenceGeneName;

    public Prediction(Match predictionResult, String name) {
        this.predictionResult = predictionResult;
        this.referenceGeneName = name;
    }

    public Match getPredictionResult() {
        return predictionResult;
    }

    public String getReferenceGeneName() {
        return referenceGeneName;
    }
}
