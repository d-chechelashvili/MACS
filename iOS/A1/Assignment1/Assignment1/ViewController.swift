//
//  ViewController.swift
//  Assignment1
//
//  Created by dato che on 10/4/20.
//

import UIKit

class ViewController: UIViewController {
    var instructions : [[Bool]] = [[true, false, true, true, true, true, true], //0
                                [false, false, false, false, true, false, true], //1
                                [true, true, true, false, true, true, false], //2
                                [true, true, true, false, true, false, true], //3
                                [false, true, false, true, true, false, true], //4
                                [true, true, true, true, false, false, true], //5
                                [true, true, true, true, false, true, true], //6
                                [false, false, true, false, true, false, true], //7
                                [true, true, true, true, true, true, true], //8
                                [true, true, true, true, true, false, true]] //9
    
    
    @IBOutlet var firstViewController: [UIView]!
    @IBOutlet var secondViewController: [UIView]!
    @IBOutlet var thirdViewController: [UIView]!
    
    var destinationNum = 342
    
    override func viewDidLoad() {
        super.viewDidLoad()
        let firstNum = destinationNum / 100
        let secondNum = (destinationNum / 10) % 10
        let thirdNum = destinationNum % 10
        
        
        viewsWillBecomeActive(num: firstNum, views: firstViewController)
        viewsWillBecomeActive(num: secondNum, views: secondViewController)
        viewsWillBecomeActive(num: thirdNum, views: thirdViewController)
    }
    
    func viewsWillBecomeActive(num:Int, views:[UIView]!){
        for i in 0...instructions[num].count - 1{
            if(instructions[num][i]){
                views[i].backgroundColor = .magenta
            }else{
                views[i].backgroundColor = .systemYellow
            }
        }
    }
}
