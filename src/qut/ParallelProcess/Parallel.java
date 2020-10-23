package qut.ParallelProcess;

import qut.*;
import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;
import java.util.concurrent.locks.ReentrantLock;
import qut.parallel.GeneComparisionTask;


public class Parallel {
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static final ThreadLocal<Series> sigma70_pattern = ThreadLocal.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];
    private static ReentrantLock lock = new ReentrantLock();


    static {
        complement['C'] = 'G';
        complement['c'] = 'g';
        complement['G'] = 'C';
        complement['g'] = 'c';
        complement['T'] = 'A';
        complement['t'] = 'a';
        complement['A'] = 'T';
        complement['a'] = 't';
    }


    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true) {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    private static boolean Homologous(PeptideSequence A, PeptideSequence B) {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location - 1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location - upStreamDistance - 1, gene.location - 1));
        else {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i = 0; i < upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart - i]];
            return new NucleotideSequence(result);
        }
    }

    private static Match PredictPromoter(NucleotideSequence upStreamRegion) {
        return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir) {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    private static List<String> ListGenbankFiles(String dir) {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    private static GenbankRecord Parse(String file) throws IOException {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }


    public void run_parallel_parse(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<GenbankRecord> records = ListGenbankFiles(dir).parallelStream()
                .map(file -> {
                    try {
                        System.out.println(file);
                        return Parse(file);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    return null;
                })
                .collect(Collectors.toList());

        for (GenbankRecord record : records) {
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    if (Homologous(gene.sequence, referenceGene.sequence)) {
                        NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                        Match prediction = PredictPromoter(upStreamRegion);
                        if (prediction != null) {
                            consensus.get(referenceGene.name).addMatch(prediction);
                            consensus.get("all").addMatch(prediction);
                        }
                    }
                }
            }
        }
    }

    public void run_parallel_homologous_parallelStream(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<GeneComparisionTask> geneComparisionTasks = new ArrayList<>();

        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    geneComparisionTasks.add(new GeneComparisionTask(record.nucleotides, referenceGene, gene));
                }
            }
        }

        List<GeneComparisionTask> comparisionResult = geneComparisionTasks.parallelStream()
                .filter(geneComparisionTask -> Homologous(geneComparisionTask.getGene().sequence, geneComparisionTask.getReferenceGene().sequence))
                .collect(Collectors.toList());

        for (GeneComparisionTask result : comparisionResult) {
            NucleotideSequence upStreamRegion = GetUpstreamRegion(result.getNucleotideSequence(), result.getGene());
            Match prediction = PredictPromoter(upStreamRegion);
            if (prediction != null) {
                consensus.get(result.getReferenceGene().name).addMatch(prediction);
                consensus.get("all").addMatch(prediction);
            }
        }
    }

    public void run_parallel_predictPromoter_parallelStream(String referenceFile, String dir) throws FileNotFoundException, IOException {

        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<GeneComparisionTask> geneComparisionTasks = new ArrayList<>();

        for (String filename : ListGenbankFiles(dir)) {
            GenbankRecord record = Parse(filename);
            System.out.println(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    geneComparisionTasks.add(new GeneComparisionTask(record.nucleotides, referenceGene, gene));
                }
            }
        }

        List<GeneComparisionTask> comparisionResult = geneComparisionTasks.parallelStream()
                .filter(geneComparisionTask -> Homologous(geneComparisionTask.getGene().sequence, geneComparisionTask.getReferenceGene().sequence))
                .map(geneComparisionTask -> {
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(geneComparisionTask.getNucleotideSequence(), geneComparisionTask.getGene());
                    Match prediction = PredictPromoter(upStreamRegion);
                    geneComparisionTask.setPrediction(prediction);
                    return geneComparisionTask;
                })
                .collect(Collectors.toList());

        for (GeneComparisionTask result : comparisionResult) {
            if (result.getPrediction() != null) {
                consensus.get(result.getReferenceGene().name).addMatch(result.getPrediction());
                consensus.get("all").addMatch(result.getPrediction());
            }
        }
    }

    public void run_parallel_addMatch_parallelStream(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<GeneComparisionTask> geneComparisionTasks = new ArrayList<>();

        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    geneComparisionTasks.add(new GeneComparisionTask(record.nucleotides, referenceGene, gene));
                }
            }
        }

        geneComparisionTasks.parallelStream()
                .filter(geneComparisionTask -> Homologous(geneComparisionTask.getGene().sequence, geneComparisionTask.getReferenceGene().sequence))
                .forEach(geneComparisionTask -> {
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(geneComparisionTask.getNucleotideSequence(), geneComparisionTask.getGene());
                    Match prediction = PredictPromoter(upStreamRegion);
                    if (prediction != null) {
                        lock.lock();
                        consensus.get(geneComparisionTask.getReferenceGene().name).addMatch(prediction);
                        consensus.get("all").addMatch(prediction);
                        lock.unlock();
                    }
                });
    }

    public void run_parallel_homologous_executorService(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        ExecutorService executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future> comparisionResult = new ArrayList<>();

        for (String filename : ListGenbankFiles(dir)) {
            GenbankRecord record = Parse(filename);
            System.out.println(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    Future<CallableTaskHomologous> futureResult = executorService.submit(new CallableTaskHomologous(record.nucleotides, referenceGene, gene));
                    comparisionResult.add(futureResult);
                }
            }
        }

        for (Future future : comparisionResult) {
            try {
                CallableTaskHomologous result = (CallableTaskHomologous) future.get();
                if (result.result) {
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(result.getNucleotideSequence(), result.getGene());
                    Match prediction = PredictPromoter(upStreamRegion);
                    if (prediction != null) {
                        consensus.get(result.getReferenceGene().name).addMatch(prediction);
                        consensus.get("all").addMatch(prediction);
                    }
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
    }

    public void run_parallel_predictPromoter_executorService(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        ExecutorService executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future> comparisionResult = new ArrayList<>();

        for (String filename : ListGenbankFiles(dir)) {
            GenbankRecord record = Parse(filename);
            System.out.println(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    Future<CallableTaskPredictPromoter> futureResult = executorService.submit(new CallableTaskPredictPromoter(record.nucleotides, referenceGene, gene));
                    comparisionResult.add(futureResult);
                }
            }
        }

        for (Future future : comparisionResult) {
            try {
                CallableTaskPredictPromoter result = (CallableTaskPredictPromoter) future.get();
                if (result.prediction != null) {
                    consensus.get(result.getReferenceGene().name).addMatch(result.prediction);
                    consensus.get("all").addMatch(result.prediction);
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
    }

    public void run_parallel_addMatch_executorService(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        ExecutorService executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        List<Future> results = new ArrayList<>();

        for (String filename : ListGenbankFiles(dir)) {
            GenbankRecord record = Parse(filename);
            System.out.println(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes) {
                    Future futureResult = executorService.submit(new RunnableTaskAddMatch(record.nucleotides, referenceGene, gene));
                    results.add(futureResult);
                }
            }
        }
        for (Future future : results) {
            try {
                future.get();
            } catch (InterruptedException e) {
                e.printStackTrace();
            } catch (ExecutionException e) {
                e.printStackTrace();
            }
        }
    }

    public void run_sequential(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        for (String filename : ListGenbankFiles(dir)) {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes) {
                System.out.println(referenceGene.name);
                for (Gene gene : record.genes){
                    if (Homologous(gene.sequence, referenceGene.sequence)) {
                        NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                        Match prediction = PredictPromoter(upStreamRegion);
                        if (prediction != null) {
                            consensus.get(referenceGene.name).addMatch(prediction);
                            consensus.get("all").addMatch(prediction);
                        }
                    }
                }
            }
        }
    }

    public class CallableTaskHomologous implements Callable<CallableTaskHomologous> {
        private final NucleotideSequence nucleotideSequence;
        private final Gene gene;
        private final Gene referenceGene;
        private Boolean result;

        public CallableTaskHomologous(NucleotideSequence nucleotideSequence, Gene referenceGene, Gene gene) {
            this.nucleotideSequence = nucleotideSequence;
            this.referenceGene = referenceGene;
            this.gene = gene;
            this.result = false;
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

        public Boolean getResult() {
            return result;
        }

        @Override
        public CallableTaskHomologous call() throws Exception {
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                this.result = true;
                return this;
            }
            return this;
        }
    }

    public class CallableTaskPredictPromoter implements Callable<CallableTaskPredictPromoter> {
        private final NucleotideSequence nucleotideSequence;
        private final Gene gene;
        private final Gene referenceGene;
        private Match prediction;

        public CallableTaskPredictPromoter(NucleotideSequence nucleotideSequence, Gene referenceGene, Gene gene) {
            this.nucleotideSequence = nucleotideSequence;
            this.referenceGene = referenceGene;
            this.gene = gene;
            this.prediction = null;
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

        public Match getPrediction() {
            return prediction;
        }

        @Override
        public CallableTaskPredictPromoter call() throws Exception {
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                NucleotideSequence upStreamRegion = GetUpstreamRegion(nucleotideSequence, gene);
                Match prediction = PredictPromoter(upStreamRegion);
                this.prediction = prediction;
                return this;
            }
            return this;
        }
    }

    public class RunnableTaskAddMatch implements Runnable {
        private final NucleotideSequence nucleotideSequence;
        private final Gene gene;
        private final Gene referenceGene;

        public RunnableTaskAddMatch(NucleotideSequence nucleotideSequence, Gene referenceGene, Gene gene) {
            this.nucleotideSequence = nucleotideSequence;
            this.referenceGene = referenceGene;
            this.gene = gene;
        }

        @Override
        public void run() {
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                NucleotideSequence upStreamRegion = GetUpstreamRegion(nucleotideSequence, gene);
                Match prediction = PredictPromoter(upStreamRegion);
                if (prediction != null) {
                    lock.lock();
                    consensus.get(referenceGene.name).addMatch(prediction);
                    consensus.get("all").addMatch(prediction);
                    lock.unlock();
                }
            }
        }
    }

    public static void main(String[] args) throws FileNotFoundException, IOException {
        long average = 0;
        String referenceFile = "referenceGenes.list"; //Modify the string with your path to the file to run the program
        String directory= "Ecoli"; //Modify the string with your path to the folder to run the program
        for (int i = 0; i < 10; i++) {
            long start = System.currentTimeMillis();
            new Parallel().run_parallel_homologous_executorService(referenceFile, directory);
            long end = System.currentTimeMillis();
            System.out.println(String.format("Run for: %s seconds", (end - start) / 1000));
            average += (end - start) / 1000;
        }
        System.out.println("Average: " + average / 10 + " seconds");
    }
}
