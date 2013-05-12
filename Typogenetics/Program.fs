// Typogenetics
// An implementation of the typogenetical system from
// Godel, Escher and Bach 504-513

type Base =
    | Adenine
    | Cytosine
    | Guanine
    | Thyamine

let strand s =
    let fromLetter = function
        | 'A' -> Adenine
        | 'C' -> Cytosine
        | 'G' -> Guanine
        | 'T' -> Thyamine
        |  _  -> failwith("Invalid Base code")
    Seq.map fromLetter s |> Seq.toList

type Codon =
    | Duplet  of Base * Base
    | Singlet of Base
    
let rec codons = function
    | x0::x1::xs -> Duplet(x0, x1)::(codons xs)
    | x::xs      -> Singlet(x)::(codons xs)
    | []         -> []
  
type AminoAcid =
    | Cut
    | Del
    | Swi
    | Mvr
    | Mvl
    | Cop
    | Off
    | Ina
    | Inc
    | Ing
    | Int
    | Rpy
    | Rpu
    | Lpy
    | Lpu
  
let decode = function
    | Duplet(Adenine,  Adenine)  -> None
    | Duplet(Adenine,  Cytosine) -> Some(Cut)
    | Duplet(Adenine,  Guanine)  -> Some(Del)
    | Duplet(Adenine,  Thyamine) -> Some(Swi)
    | Duplet(Cytosine, Adenine)  -> Some(Mvr)
    | Duplet(Cytosine, Cytosine) -> Some(Mvl)
    | Duplet(Cytosine, Guanine)  -> Some(Cop)
    | Duplet(Cytosine, Thyamine) -> Some(Off)
    | Duplet(Guanine,  Adenine)  -> Some(Ina)
    | Duplet(Guanine,  Cytosine) -> Some(Inc)
    | Duplet(Guanine,  Guanine)  -> Some(Ing)
    | Duplet(Guanine,  Thyamine) -> Some(Int)
    | Duplet(Thyamine, Adenine)  -> Some(Rpy)
    | Duplet(Thyamine, Cytosine) -> Some(Rpu)
    | Duplet(Thyamine, Guanine)  -> Some(Lpy)
    | Duplet(Thyamine, Thyamine) -> Some(Lpu)
    | Singlet(_)                 -> None


let translate =
    List.map decode

let split delimiter list =
    let rec split' list' =
        match list' with
        | []                       -> [[]]
        | x::xs when x = delimiter -> [] :: split' xs
        | x::xs                    -> let sxs = split' xs
                                      (x :: sxs.Head) :: sxs.Tail
    split' list

let splitn n ls =
    let rec split' i source prefix suffix =
        match source with
        | []               -> (prefix, suffix)
        | x::xs when i < n -> let p, s = split' (i + 1) xs prefix suffix
                              (x::p, s)      
        | x::xs            -> let p, s = split' (i + 1) xs prefix suffix
                              (p, x::s)
    split' 0 ls [] [] 

type Direction =
    | Straight
    | Left
    | Right

let direction = function
    | Cut -> Straight
    | Del -> Straight
    | Swi -> Right
    | Mvr -> Straight
    | Mvl -> Straight
    | Cop -> Right
    | Off -> Left
    | Ina -> Straight
    | Inc -> Right
    | Ing -> Right
    | Int -> Left
    | Rpy -> Right
    | Rpu -> Left
    | Lpy -> Left
    | Lpu -> Left
    
let s = "GATCGAATGTAGAGTGATGCTGCTGCTAGTACTGA"

let stripOptions list =
    List.map Option.get list

let turn = function
    | Left     ->  90
    | Straight ->   0
    | Right    -> -90
        
let intern (enzyme:'a list) =
    if enzyme.Length < 3 then
        Seq.empty
    else
        Seq.take (enzyme.Length - 2) (enzyme.Tail)

let bindingPreference enzyme =
    let truncated = intern enzyme
    let angle = Seq.map direction truncated |> Seq.map turn |> Seq.fold (+) 0
    match angle % 360 with
    |   0 -> Adenine
    |  90 -> Cytosine
    | 180 -> Guanine
    | 270 -> Thyamine
    |  _  -> failwith("Unexpected angle") 

let dna = s|> strand

let enzymes = dna |> codons |> translate |> split None |> List.filter (fun x -> not x.IsEmpty) |> List.map stripOptions

type Mode =
    | Move
    | Copy
    | Quit

let complementBase b =
    match b with
    | Adenine  -> Thyamine
    | Cytosine -> Guanine
    | Guanine  -> Cytosine
    | Thyamine -> Adenine

// If mode is true, affirms that working[index]
// is copied (in the complementary pairing sense) to
// complement[index]
let affirmCopy working complement index mode =
    let complement' =
        match mode with
        | Quit -> complement
        | Move -> complement
        | Copy -> let prefix, suffix = splitn index complement
                  let b = Option.map complementBase (List.nth working index)
                  prefix @ [ b ] @ suffix.Tail    
    (working, complement', index, mode)

// Cut both strands to the right of the index position. Any fragments to the
// right are added to leftovers.
let cutStrands working complement leftovers index mode =
    let cutStrand strand =
        let keep, remainder = splitn (index + 1) strand
        let leftovers = split None remainder
        (keep, leftovers)
    let working', workingLeftovers = cutStrand working
    let complement', complementLeftovers = cutStrand complement
    let leftovers' = workingLeftovers @ complementLeftovers
    
    (working', complement', leftovers', index, Quit)
    
// Removes the base at the current index from working,
// leaving its complement untouched
let deleteBaseFromStrand working complement leftovers index mode =
    let working' = List.mapi (fun i b -> if i = index then None else b) working
    (working', complement, leftovers, index + 1, mode) 

let switchToOtherStrand working complement leftovers index mode =
    let working'    = List.rev complement
    let complement' = List.rev working
    (working', complement', leftovers, index, mode)

let move direction working complement leftovers index mode =
    let working', complement', index', mode' =
            affirmCopy working complement index mode
    (working', complement', leftovers, (direction index), mode)

let modeOn working complement leftovers index mode =
    let working', complement', index', mode' =
        affirmCopy working complement index Copy
    (working', complement', leftovers, index, mode)

let modeOff working complement leftovers index mode =
    (working, complement, leftovers, index, Move)

let insert b working complement leftovers index mode  =
    let workingPrefix, workingSuffix = splitn index working
    let newWorkingBase = Some(b)
    let working' = workingPrefix @ [ newWorkingBase ] @ workingSuffix
    let working', complement', index', mode' =
        affirmCopy working' complement index mode
    (working', complement', leftovers, index, mode)
    
let (|Purine|_|) b =
    match b with
    | Adenine | Guanine -> Some()
    | _                 -> None
    
let (|Pyrimidine|_|) b =
    match b with
    | Cytosine | Thyamine -> Some()
    | _                   -> None
    
let right index = index - 1
let left index  = index + 1
    
let search (|Pat|_|) direction working complement leftovers index mode =
    let rec search' working complement index mode =
        let working', complement', index', mode' =
            affirmCopy working complement index mode
        match List.nth working' index' with
        | Some(Pat)        -> (working', complement', index', mode')
        | Some(_)          -> search' working' complement' (direction index') mode'
        | None             -> (working', complement', index', Quit)
    let working'', complement'', index'', mode'' =
        search' working complement index mode
    (working'', complement'', leftovers, index'', mode'')

// Apply a single amino acid to the strands at the specified index
// and returns the modifies strands
let applyAcid (working, complement, leftovers, index, mode) acid =
    let applicator =
        match acid with
        | Cut -> cutStrands
        | Del -> deleteBaseFromStrand
        | Swi -> switchToOtherStrand
        | Mvr -> move right
        | Mvl -> move left
        | Cop -> modeOn
        | Off -> modeOff
        | Ina -> insert Adenine
        | Inc -> insert Cytosine
        | Ing -> insert Guanine
        | Int -> insert Thyamine
        | Rpy -> search (|Pyrimidine|_|) right
        | Rpu -> search (|Purine|_|)     right
        | Lpy -> search (|Pyrimidine|_|) left
        | Lpu -> search (|Purine|_|)     left
    applicator working complement leftovers index mode

let unzipStrands strands =
    let splitDoubleStrands (working, complement, leftovers) =
        [working; List.rev complement] @ leftovers
    let splitStrand strand =
        split None strand  
    strands |> splitDoubleStrands |> List.map splitStrand |> List.concat |> List.filter (fun x -> not x.IsEmpty)
 
// Given an enzyme, apply it to a strand, with a ribosome
let applyEnzyme enzyme strand =
    // Copy the strand into an array
    let workingStrand = List.map Some strand
    
    // Create an array for the complementary strand
    let complementaryStrand = List.replicate (List.length workingStrand) Option<Base>.None
       
    // Bind the enzyme to the strand - set the index
    let basePredicate unit =
        match unit with
        | Some(nucleobase) -> nucleobase = bindingPreference enzyme
        | None             -> false
    
    let startIndex = List.tryFindIndex basePredicate workingStrand
      
    // Process each amino acid from the enzyme
    let synthesize working complement index =
        let rec synthesize' (w, c, l, index, mode) enz =
            printf "%A\n" ((w, c, l, index, mode), enz)
            match mode, enz with
            | Quit, _           -> (w, c, l)
            | _   , acid::acids -> synthesize' (applyAcid (w, c, l, index, mode) acid) acids
            | _   , []          -> (w, c, l)
        synthesize' (working, complement, [], index, Move) enzyme     
                   
    let apply =
        match startIndex with
        | Some(i) -> synthesize workingStrand complementaryStrand i
        | None    -> (workingStrand, complementaryStrand, [])

    apply |> unzipStrands

[<EntryPoint>]
let main argv = 
    printfn "%A" argv
    0 // return an integer exit code
