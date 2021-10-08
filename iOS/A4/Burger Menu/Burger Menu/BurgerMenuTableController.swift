//
//  BurgerMenuTableController.swift
//  Burger Menu
//
//  Created by dato che on 12/8/20.
//

import UIKit


class BurgerMenuSection{
    var header: BurgerMenuHeaderModel?
    var burgerMenuButtons = [ButtonCellModel]()
    
    var numberOfRows: Int{
        if let header = header{
            return header.isExpanded ? burgerMenuButtons.count : 0
        }
        return burgerMenuButtons.count
    }
    
    init(header: BurgerMenuHeaderModel?, burgerMenuButtons: [ButtonCellModel]) {
        self.header = header
        self.burgerMenuButtons = burgerMenuButtons
    }
}

class BurgerMenuTableController: UIViewController {
    
    @IBOutlet var burgerMenuTableView: UITableView!

    private var tableData = [BurgerMenuSection]()
    
    override func viewDidLoad() {
        super.viewDidLoad()
        
        initData()
        
        burgerMenuTableView.delegate = self
        burgerMenuTableView.dataSource = self
        
        burgerMenuTableView.register(
            UINib(nibName: "ButtonCell", bundle: nil),
            forCellReuseIdentifier: "ButtonCell")
        
        burgerMenuTableView.register(
            UINib(nibName: "BurgerMenuHeader", bundle: nil),
            forHeaderFooterViewReuseIdentifier: "BurgerMenuHeader")

        burgerMenuTableView.tableFooterView = UIView(frame: .zero)
        burgerMenuTableView.backgroundColor = .cyan
    }
    
    private func initData(){
        tableData = [
            BurgerMenuSection(
                header: BurgerMenuHeaderModel(isExpanded: false, delegate: self),
                burgerMenuButtons: [ ButtonCellModel(text: "about us"),
                                     ButtonCellModel(text: "products"),
                                     ButtonCellModel(text: "media"),
                                     ButtonCellModel(text: "contact us")]
            )
        ]

    }
}


extension BurgerMenuTableController: UITableViewDelegate, UITableViewDataSource{
    
    func numberOfSections(in tableView: UITableView) -> Int {
        return tableData.count
    }
    
    func tableView(_ tableView: UITableView, numberOfRowsInSection section: Int) -> Int {
        return tableData[section].numberOfRows
    }
    
    func tableView(_ tableView: UITableView, cellForRowAt indexPath: IndexPath) -> UITableViewCell {
        let cell = burgerMenuTableView.dequeueReusableCell(withIdentifier: "ButtonCell", for: indexPath)
        if let buttonCell = cell as? ButtonCell{
            buttonCell.configure(with: tableData[indexPath.section].burgerMenuButtons[indexPath.row])
        }
        return cell
    }
    
    func tableView(_ tableView: UITableView, didSelectRowAt indexPath: IndexPath) {
        burgerMenuTableView.deselectRow(at: indexPath, animated: true)
    }
    
    func tableView(_ tableView: UITableView, viewForHeaderInSection section: Int) -> UIView? {
        guard let headerModel = tableData[section].header else {return nil}
        let header = burgerMenuTableView.dequeueReusableHeaderFooterView(withIdentifier: "BurgerMenuHeader")
        // 93-e xazze wers aset rames - "Changing the background color of UITableViewHeaderFooterView is not supported. Use the background view configuration instead."
        //warmodgena ar maqvs risi bralia, me ar vcvli arsad UITableViewHeaderFooterView -s background colors
        if let burgerMenuHeader = header as? BurgerMenuHeader{
            burgerMenuHeader.configure(with: headerModel)
        }
        return header
    }
    
    func tableView(_ tableView: UITableView, heightForHeaderInSection section: Int) -> CGFloat {
        return burgerMenuTableView.frame.height * 0.16
    }
    
    func tableView(_ tableView: UITableView, heightForRowAt indexPath: IndexPath) -> CGFloat {
        return burgerMenuTableView.frame.height * 0.085
    }
}


extension BurgerMenuTableController: BurgerMenuHeaderDelegate{
    func BurgerMenuHeaderDidClickButton(_ header: BurgerMenuHeader) {
        burgerMenuTableView.beginUpdates()
        if burgerMenuTableView.backgroundColor == UIColor.darkGray{
            burgerMenuTableView.backgroundColor = .cyan
        }else{
            burgerMenuTableView.backgroundColor = .darkGray
        }
        burgerMenuTableView.reloadSections(IndexSet(integer: 0), with: .top)
        burgerMenuTableView.endUpdates()
    }
}
