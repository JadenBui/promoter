package qut.parallel;

import qut.*;
import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import qut.parallel.EqualityTest;

import java.io.*;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;


public class Parallel
{
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static final ThreadLocal<Series> sigma70_pattern = ThreadLocal.withInitial(() ->Sigma70Definition.getSeriesAll_Unanchored(0.7));
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];
    private static ReentrantLock lock = new ReentrantLock();
    public static String version = null;

    static
    {
        complement['C'] = 'G'; complement['c'] = 'g';
        complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a';
        complement['A'] = 'T'; complement['a'] = 't';
    }

    public static HashMap<String, Sigma70Consensus> getConsensus() {
        return consensus;
    }

    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true)
        {
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

    private static boolean Homologous(PeptideSequence A, PeptideSequence B)
    {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene)
    {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location-1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
        else
        {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i=0; i<upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart-i]];
            return new NucleotideSequence(result);
        }
    }

    private static Match PredictPromoter(NucleotideSequence upStreamRegion)
    {
        return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir)
    {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    private static List<String> ListGenbankFiles(String dir)
    {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    private static GenbankRecord Parse(String file) throws IOException
    {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    public void run_parallel_executorServices(String referenceFile, String dir) throws FileNotFoundException, IOException
    {
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

    public void run_parallel_parallelStream(String referenceFile, String dir) throws FileNotFoundException, IOException {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        List<GeneComparisionTask> geneComparisionTasks = new ArrayList<>();
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "12");

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
        public void run(){
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                NucleotideSequence upStreamRegion = GetUpstreamRegion(nucleotideSequence, gene);
                Match prediction = PredictPromoter(upStreamRegion);
                if (prediction != null)
                {
                    lock.lock();
                    consensus.get(referenceGene.name).addMatch(prediction);
                    consensus.get("all").addMatch(prediction);
                    lock.unlock();
                }
            }
        }
    }

    public static void main(String[] args) throws FileNotFoundException, IOException
    {
        Scanner scanner = new Scanner(System.in);
        while(version == null){
            System.out.print("Please enter the parallel version: \n"+
                    "1. Executors API \n"+
                    "2. Parallel Stream API \n" +
                    "Your choice (1|2): ");
            int choice = scanner.nextInt();
            version = choice == 1 ? "Executors API" : "Parallel Stream API";
        }
        long start = System.currentTimeMillis();
        if(version.equals("Executors API")){
            new Parallel().run_parallel_executorServices("C:\\Users\\sheep\\CAB401_Program\\promoter\\referenceGenes.list", "C:\\Users\\sheep\\CAB401_Program\\promoter\\Ecoli");
        }else{
            new Parallel().run_parallel_parallelStream("C:\\Users\\sheep\\CAB401_Program\\promoter\\referenceGenes.list", "C:\\Users\\sheep\\CAB401_Program\\promoter\\Ecoli");
        }
        long end = System.currentTimeMillis();
        System.out.println(String.format(version+" version ran for: %s seconds", (end - start) / 1000));
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
        EqualityTest.assertEquals(consensus);
    }
}
