//
//  BurgerMenuHeader.swift
//  Burger Menu
//
//  Created by dato che on 12/8/20.
//

import UIKit

class BurgerMenuHeaderModel{
    var isExpanded: Bool
    weak var delegate: BurgerMenuHeaderDelegate?

    init(isExpanded: Bool, delegate: BurgerMenuHeaderDelegate) {
        self.isExpanded = isExpanded
        self.delegate = delegate
    }
}

protocol BurgerMenuHeaderDelegate: AnyObject {
    func BurgerMenuHeaderDidClickButton(_ header: BurgerMenuHeader)
}

class BurgerMenuHeader: UITableViewHeaderFooterView{
    
    @IBOutlet var menuButton: UIButton!
    
    var model: BurgerMenuHeaderModel!
    
    func configure(with model: BurgerMenuHeaderModel){
       //virtualboxis gachedvis gamo ver vigeb animacia ramdenad sworad mushaobs
        //magram memgoni sworad mushaobs :D
        UIView.animate(
            withDuration: 1,
            animations: {
                if model.isExpanded{
                    self.menuButton.transform = CGAffineTransform.init(rotationAngle: -.pi)
                    self.menuButton.setImage(UIImage(systemName: "x.circle"), for: .normal)
                }else{
                    self.menuButton.transform = CGAffineTransform.init(rotationAngle: .pi)
                    self.menuButton.setImage(UIImage(systemName: "line.horizontal.3"), for: .normal)
                }
            })
        self.model = model
    }
    
    @IBAction func handleMenuButton(){
        model.isExpanded.toggle()
        model.delegate?.BurgerMenuHeaderDidClickButton(self)
    }
}
