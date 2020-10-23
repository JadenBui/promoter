package qut.parallel;

import org.junit.jupiter.api.Test;

import java.io.IOException;

import static org.junit.jupiter.api.Assertions.*;

class ParallelTest {
    /*Test Parallel Stream API version*/
    @Test
    void parallelStreamAPI() {
        try {
            Parallel.version = "Parallel Stream API";
            Parallel.main(new String[0]);
            assertTrue(EqualityTest.assertEqualsBoolean(Parallel.getConsensus()));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    /*Test Executors API version*/
    @Test
    void ExecutorsAPI() {
        try {
            Parallel.version = "Executors API";
            Parallel.main(new String[0]);
            assertTrue(EqualityTest.assertEqualsBoolean(Parallel.getConsensus()));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}