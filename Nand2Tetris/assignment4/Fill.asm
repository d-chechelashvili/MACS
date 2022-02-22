// This file is part of www.nand2tetris.org
// and the book "The Elements of Computing Systems"
// by Nisan and Schocken, MIT Press.
// File name: projects/04/Fill.asm

// Runs an infinite loop that listens to the keyboard input.
// When a key is pressed (any key), the program blackens the screen,
// i.e. writes "black" in every pixel;
// the screen should remain fully black as long as the key is pressed. 
// When no key is pressed, the program clears the screen, i.e. writes
// "white" in every pixel;
// the screen should remain fully clear as long as no key is pressed.

// Put your code here.

// while(true){
//     addr = SCREEN
//     n = 8192
//     i = 0
//     while (i < n){
//         if(KBD){
//             RAM[SCREEN + 16*i] = -1
//         }else{
//             RAM[SCREEN + 16*i] = 0
//         }
//         addr++
//         i++
//     } 
// }


(WHILETRUE)
    @SCREEN
    D=A
    @addr
    M=D     // addr = 16384 (screen's base address)

    @8192
    D=A
    @n
    M=D     // n = 8192 (screen length/16)

    @i
    M=0     // i = 0
(LOOP)
    //i < n
    @i
    D=M
    @n
    D=D-M
    @WHILETRUE
    D;JEQ

    @KBD
    D=M

    @NOTPRESSED 
    D;JEQ   // if keyboard is not pressed jump

    @addr
    A=M
    M=-1    // set pixels to 1
    @ENDIF
    0;JMP   //goto end of if-else
(NOTPRESSED)
    @addr
    A=M
    M=0     // set pixels to 0 
(ENDIF)
    @addr
    M=M+1  //addr++
    @i
    M=M+1   //i++
    @LOOP
    0;JMP   //goto loop


