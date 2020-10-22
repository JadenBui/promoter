package qut.parallel;

import com.google.gson.Gson;
import qut.Sigma70Consensus;

import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;

public class EqualityTest {
    public static void assertEquals(Map<String, Sigma70Consensus> consensus) {
        try {
            Gson gson = new Gson();
            Reader reader = Files.newBufferedReader(Paths.get("consensus.json"));
            String resultInJson = gson.toJson(consensus);
            Map<String, Sigma70Consensus> originalResult = gson.fromJson(reader, Map.class);
            String originalResultToString = gson.toJson(originalResult).replace(".0","");
            System.out.println("The result after parallelism is " +
                    (resultInJson.equals(originalResultToString) ? "equal" : "not equal")
            + " to the original result.");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static boolean assertEqualsBoolean(Map<String, Sigma70Consensus> consensus) {
        try {
            Gson gson = new Gson();
            Reader reader = Files.newBufferedReader(Paths.get("consensus.json"));
            String resultInJson = gson.toJson(consensus);
            Map<String, Sigma70Consensus> originalResult = gson.fromJson(reader, Map.class);
            String originalResultToString = gson.toJson(originalResult).replace(".0","");
            return resultInJson.equals(originalResultToString);
        } catch (Exception ex) {
            ex.printStackTrace();
            throw new Error();
        }
    }
}
