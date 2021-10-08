//
//  ViewController.swift
//  Bidirectional Layout
//
//  Created by dato che on 12/15/20.
//

import UIKit

class ViewController: UIViewController {

    @IBOutlet var BidirectionalTableView: UITableView!

    override func viewDidLoad() {
        super.viewDidLoad()
        
        BidirectionalTableView.delegate = self
        BidirectionalTableView.dataSource = self
        
        BidirectionalTableView.separatorStyle = .none
        BidirectionalTableView.allowsSelection = false
        
        BidirectionalTableView.register(UINib(nibName: "FirstSectionCell", bundle: nil),
                               forCellReuseIdentifier: "FirstSectionCell")
        BidirectionalTableView.register(UINib(nibName: "FullRedCell", bundle: nil),
                               forCellReuseIdentifier: "FullRedCell")
        BidirectionalTableView.register(UINib(nibName: "CollectionViewAsCell", bundle: nil),
                               forCellReuseIdentifier: "CollectionViewAsCell")
        
        BidirectionalTableView.backgroundColor = .white
        
    }
}

extension ViewController: UITableViewDelegate, UITableViewDataSource{
    
    func numberOfSections(in tableView: UITableView) -> Int {
        return 7
    }
    
    func tableView(_ tableView: UITableView, numberOfRowsInSection section: Int) -> Int {
        return 1
    }
    
    func tableView(_ tableView: UITableView, cellForRowAt indexPath: IndexPath) -> UITableViewCell {
        if indexPath.section == 0{
            let cell = BidirectionalTableView.dequeueReusableCell(withIdentifier: "FirstSectionCell", for: indexPath)
            return cell
        }
        if indexPath.section == 2{
            let cell = BidirectionalTableView.dequeueReusableCell(withIdentifier: "CollectionViewAsCell",
                                                                  for:indexPath)
            return cell
        }
        let cell = BidirectionalTableView.dequeueReusableCell(withIdentifier: "FullRedCell", for: indexPath)
        return cell
        
    }
    
    func tableView(_ tableView: UITableView, heightForRowAt indexPath: IndexPath) -> CGFloat {
        if indexPath.section == 0{
            return BidirectionalTableView.frame.height * 0.2
        }
        if indexPath.section == 2{
            return BidirectionalTableView.frame.height * 0.25
        }
        return BidirectionalTableView.frame.height * 0.123
    }
}
