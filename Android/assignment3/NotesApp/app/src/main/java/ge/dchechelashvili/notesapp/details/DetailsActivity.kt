package ge.dchechelashvili.notesapp.details

import android.content.Context
import android.content.Intent
import androidx.appcompat.app.AppCompatActivity
import android.os.Bundle
import android.util.Log
import android.view.View
import android.widget.EditText
import android.widget.ImageView
import androidx.recyclerview.widget.RecyclerView
import ge.dchechelashvili.notesapp.R
import ge.dchechelashvili.notesapp.data.entity.Item
import ge.dchechelashvili.notesapp.data.entity.Note

class DetailsActivity : AppCompatActivity(), IDetailsView, CheckListListener {

    private var isPinned = false
    private lateinit var presenter: DetailsPresenter
    private lateinit var rvChecked: RecyclerView
    private lateinit var rvUnchecked: RecyclerView
    private var checkedAdapter = CheckListAdapter(this)
    private var uncheckedAdapter = CheckListAdapter(this)
    private lateinit var noteName: EditText
    private var checkedItems = mutableListOf<Item>()
    private var uncheckedItems = mutableListOf<Item>()
    private var noteId = -1


    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        setContentView(R.layout.activity_details)
        presenter = DetailsPresenter(this)
        initView()
        val id = intent.getIntExtra(NOTE_ID_EXTRA, 0)
        if (id != -1){
            presenter.getNote(id)
        }
    }

    override fun onDestroy() {
        super.onDestroy()
    }

    private fun initView() {
        rvChecked = findViewById(R.id.rvChecked)
        rvChecked.adapter = checkedAdapter
        rvUnchecked = findViewById(R.id.rvUnchecked)
        rvUnchecked.adapter = uncheckedAdapter
        noteName = findViewById(R.id.editTextNoteName)
    }

    fun togglePin(view: View) {
        val pinView = view as ImageView

        if (isPinned){
            pinView.setImageResource(R.drawable.ic_pin)
            isPinned = false
        }else{
            pinView.setImageResource(R.drawable.ic_pinned)
            isPinned = true
        }
    }

    override fun onBackPressed() {
        goBack(null)
    }

    fun goBack(view: View?){
        presenter.saveInfo(noteId, noteName.text.toString(), isPinned, checkedItems, uncheckedItems)
    }

    companion object{
        const val NOTE_ID_EXTRA = "noteId"
    }

    override fun onCheckBoxClicked(item: Item, pos: Int) {
        if(item.checked){
            checkedItems.removeAt(pos)
            val newItem = Item(item.content, item.noteId, false)
            uncheckedItems.add(newItem)

        }else{
            uncheckedItems.removeAt(pos)
            val newItem = Item(item.content, item.noteId, true)
            checkedItems.add(newItem)
        }
        checkedAdapter.items = checkedItems
        uncheckedAdapter.items = uncheckedItems
        uncheckedAdapter.notifyDataSetChanged()
        checkedAdapter.notifyDataSetChanged()
    }

    override fun onRemoveButtonClicked(item: Item, pos: Int) {
        if(item.checked){
            checkedItems.removeAt(pos)
        }else{
            uncheckedItems.removeAt(pos)
        }
        checkedAdapter.items = checkedItems
        uncheckedAdapter.items = uncheckedItems
        uncheckedAdapter.notifyDataSetChanged()
        checkedAdapter.notifyDataSetChanged()
    }

    fun addItem(view: View) {
        val item = Item("",-1, false)
        uncheckedItems.add(item)
        uncheckedAdapter.items = uncheckedItems
        uncheckedAdapter.notifyDataSetChanged()
    }

    override fun exit() {
        finish()
    }

    override fun showNote(note: Note) {
        noteName.setText(note.name)
        isPinned = !note.pinned
        noteId = note.id
        togglePin(findViewById(R.id.pinView))
    }

    override fun showCheckedItems(checkedItems: MutableList<Item>) {
        this.checkedItems = checkedItems
        checkedAdapter.items = checkedItems
        checkedAdapter.notifyDataSetChanged()
    }

    override fun showUncheckedItems(uncheckedItems: MutableList<Item>) {
        this.uncheckedItems = uncheckedItems
        uncheckedAdapter.items = uncheckedItems
        uncheckedAdapter.notifyDataSetChanged()
    }

}

interface CheckListListener {
    fun onCheckBoxClicked(item: Item, pos: Int)
    abstract fun onRemoveButtonClicked(item: Item, pos: Int)
}